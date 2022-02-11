#include "Galaxy.h"
double pi = 3.141592654;
Galaxy::Galaxy(InitialisedData & data): Data(data), Param(data.Param), IGM(GasReservoir::Primordial(data.Param.Galaxy.IGM_Mass,data.Param))
{
	int currentRings = 0;
	
	Data.UrgentLog("\tMain Galaxy Initialised.\n\tStarting ring population:  ");
	
	double initMass = 0;
	double initialScaleLength = GasScaleLength(0);
	RingMasses.resize(Param.Galaxy.RingCount);
	for (int i = 0; i < Param.Galaxy.RingCount; ++i)
	{
		double ri = Param.Galaxy.RingRadius[i];
		double ringWidth = Param.Galaxy.RingWidth[i];
		double predictedDensity = PredictSurfaceDensity(ri,ringWidth,Param.Galaxy.PrimordialMass,initialScaleLength);
		double predictedMass = 2*pi * ri * ringWidth * predictedDensity;
		Rings.push_back(Ring(i,predictedMass,Data));
		Data.ProgressBar(currentRings, i,Param.Galaxy.RingCount);
		initMass += predictedMass;
	}
	
	Threads.resize(Param.Meta.ParallelThreads-1);
	Data.UrgentLog("\tGalaxy Rings initialised.\n");
	
	
	Migrator = std::vector<MigrationMatrix>(Param.Meta.SimulationSteps,MigrationMatrix(Data));
	Data.UrgentLog("\tMigration Matrices initialised.\n");
}

void Galaxy::RingEvolve(int timestep,int ringStart, int ringEnd)
{
	for (int i = ringStart; i < ringEnd; ++i)
	{
		Rings[i].TimeStep(timestep);
	}
	
}

void Galaxy::LaunchParallelOperation(int timestep, int nOperations, ParallelJob type)
{
	//~ std::cout << "A new parallel job started at " << timestep << " jobtype = " << type << "  calling " << nOperations << std::endl;
	int N = Param.Meta.ParallelThreads;
	int chunkDivisor = ceil((double)nOperations / N);
	int start = 0;
	int end = 0;
	
	int n = 0;

	while (n < N-1 && end < nOperations)
	{
		start = n*chunkDivisor;
		end = std::min(start + chunkDivisor,nOperations);
		//~ std::cout << "\tSending " << start << "->" << end << " to worker " << n << std::endl;
		switch (type)
		{
			case RingStep:
			{
				Threads[n] = std::thread(&Galaxy::RingEvolve,this,timestep,start,end);
				break;
			}
			case Scattering:
			{
				Threads[n] = std::thread(&Galaxy::ScatterYields,this,timestep,start,end);
				break;
			}
			case Compounding:
			{
				Threads[n] = std::thread(&Galaxy::CompoundScattering,this,timestep,start,end);
				break;
			}
		}
		++n;
		
	}
	
	start = end;
	end = nOperations;
	//~ std::cout << "I am working on the final set : " << start << "->" << end << std::endl;
	switch (type)
	{
		case RingStep:
		{
			RingEvolve(timestep, start,end);
			break;
		}
		case Scattering:
		{
			ScatterYields(timestep,start,end);
			break;
		}
		case Compounding:
		{
			CompoundScattering(timestep,start,end);
			break;
		}
		
	}
	
	int joined = 0;
	
	while (joined < n)
	{	
		for (int nn = 0; nn < n; ++nn)
			{
				if (Threads[nn].joinable())
				{
					Threads[nn].join();
					++joined;
				}
			}
	}
}
void Galaxy::Evolve()
{
	double t = 0;
	SaveState(t);

	int fullBar = Param.Meta.ProgressHashes;
	int currentBars = 0;
	
	Data.UrgentLog("\tStarting Galaxy evolution: ");
	
	for (int timestep = 0; timestep < Param.Meta.SimulationSteps-1; ++timestep)
	{
		//~ std::cout << "Time " << timestep << std::endl;
		IGM.PassiveCool(Param.Meta.TimeStep,true);
		//~ std::cout << "\tPassive cooled" <<std::endl;
		Infall(t);
		//~ std::cout << "\tInfell" << "  " << Mass() << std::endl;
		ComputeScattering(timestep);
		//~ std::cout << "\tComputed scattering" << "  " << Mass() << std::endl;
		LaunchParallelOperation(timestep,Rings.size(),RingStep);
		//~ std::cout << "\tRingstepped" << "  " << Mass() << std::endl;
		
		LaunchParallelOperation(timestep,Rings.size(),Scattering);
		//~ std::cout << "\tScattered yields" << "  " << Mass() << std::endl;
		ScatterGas(timestep);
		//~ std::cout << "\tScattered Gas" << "  " << Mass() << std::endl;
		t += Param.Meta.TimeStep;
		
		
		Data.ProgressBar(currentBars, timestep,Param.Meta.SimulationSteps);	
		SaveState(t);
		//~ std::cout << "\tSaved" <<std::endl;	
	}
	
	
	
}
double Galaxy::InfallMass(double t)
{
	double delta = Param.Meta.TimeStep;
	//I have analytically integrated Mdot between t-Delta and t, just for an additional layer of accuracy at early times
	double bFast = Param.Galaxy.InfallTime1;
	double bSlow = Param.Galaxy.InfallTime2;
	double fastInfall = Param.Galaxy.InfallMass1 * exp(-t/bFast) * (exp(delta / bFast) - 1.0);
	double slowInfall = Param.Galaxy.InfallMass2 * exp(-t/bSlow) * ( exp(delta/bSlow) - 1.0);
	
	
	return fastInfall + slowInfall;
}
double Galaxy::GasMass()
{
	double m = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		m+=Rings[i].Gas.Mass();
	}
	return m;
}
double Galaxy::ColdGasMass()
{
	double m = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		m+=Rings[i].Gas.ColdMass();
	}
	return m;
}
double Galaxy::StarMass()
{
	double m = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		m += Rings[i].Stars.AliveMass();
	}
	return m;
}
double Galaxy::RelicMass()
{
	double m = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		MassReport Mrr = Rings[i].Stars.DeadMass();
		m += Mrr.Total/1e9;
	}
	return m;
}
double Galaxy::Mass()
{
	return GasMass() + StarMass() + RelicMass();
}
double Galaxy::GasScaleLength(double t)
{
	
	double t0 = Param.Galaxy.ScaleLengthDelay;
	double tg = Param.Galaxy.ScaleLengthTimeScale;
	double tf = Param.Galaxy.ScaleLengthFinalTime;
	double R0 = Param.Galaxy.MinScaleLength;
	double Rf = Param.Galaxy.MaxScaleLength;

	double N = 1.0/(atan((tf - t0)/tg) - atan(-t0/tg));
	return R0 + N*(Rf - R0) * (atan( (t - t0)/tg) - atan(-t0/tg));
}

double bilitewskiRatio(double a, double b, double radius, double width, double nextwidth,double maxRadius)
{
	//~ double denominator = 1.0/(2*n + 1);
	//~ double nsq = n*n;
	//~ double firstTerm = (-2 * a)/(4 * N) * denominator * (4*n*nsq + 6*nsq + 4*n + 1);
	//~ double secondTerm = 2*(1 - b)/3 * denominator * (3 *nsq + 3*n + 1);
	
	//~ return firstTerm + secondTerm;
	
	double inflowFactor = (width + nextwidth)/2;
	double rPlus = radius + width/2;
	double rMinus = radius - width/2;
	double onflowFactor = 1.0/(radius * width);
	
	double quarticTerm = a/(4 * maxRadius) * (pow(rPlus,4) - pow(rMinus,4));
	double cubicTerm = (b - 1)/3 * (pow(rPlus,3) - pow(rMinus,3));
	
	double onflowTerm = onflowFactor * (quarticTerm + cubicTerm);
	
	return -onflowTerm / inflowFactor;
	
}

void Galaxy::InsertInfallingGas(int ring, double amount)
{
	double oldMass = Rings[ring].Gas.Mass();
	double a_factor  = Param.Migration.InflowParameterA;
	double b_factor = Param.Migration.InflowParameterB;
	double remainingMass;
	double ratio = 1;
	if ( ring < Rings.size() - 1)
	{
		double radius = Rings[ring].Radius;
		double width = Rings[ring].Width;
		double nextwidth = Rings[ring+1].Width;
		ratio = bilitewskiRatio(a_factor,b_factor,radius,width,nextwidth,Param.Galaxy.Radius);
		double inflowMass = ratio/(1 + ratio) * amount;
		//check that we do not remove more gas than is actually present
		double maxDepletion = Param.Migration.MaxStealFraction;
		inflowMass = std::min(inflowMass, maxDepletion*Rings[ring+1].Gas.Mass());
		//~ std::cout << a_factor << "  " << b_factor << "  " << inflowMass <<std::endl;
				
		Rings[ring].Gas.TransferFrom(Rings[ring+1].Gas,inflowMass);
		//if some part of the budget was missed because of the std::min above, then make up the deficit from the IGM
		remainingMass = amount -inflowMass;
	}
	else
	{
		remainingMass = amount;
	}
	Rings[ring].Gas.Absorb(IGM.AccretionStream(remainingMass));
}

std::vector<double> IterativeFit(const std::vector<double> & oldDeltas, const double newMass)
{
	//~ std::cout << "Had to use iterative" <<std::endl;
	int n = oldDeltas.size();
	std::vector<double> newDeltas(n,0.0);
	double correctedAmount = 0;
	for (int i = 0; i < n; ++i)
	{
		double proposal = std::max(0.0,oldDeltas[i]);
		correctedAmount += proposal;
		newDeltas[i] = proposal;
	}
	
	//~ double correctionFactor = 0;
	//~ if (newMass > 0)
	//~ {
	double correctionFactor = newMass / (correctedAmount);
	//~ }
	double sumsum = 0;
	if (!std::isnan(correctionFactor))
	{
		for (int i = 0; i < n; ++i)
		{
			newDeltas[i] *= correctionFactor;
			sumsum += newDeltas[i];
		}
	}
	return newDeltas;
	
}

void Galaxy::Infall(double t)
{
	double Rd = GasScaleLength(t);
	double oldGas = GasMass();
	double predictedInfall = InfallMass(t);
	
	double newGas = oldGas + predictedInfall;
	
	std::vector<double> origMass(Rings.size(),0.0);
	std::vector<double> perfectMasses(Rings.size(),0.0);
	std::vector<double> perfectDeltas(Rings.size(),0.0);
	bool perfect = true;
	double perf = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		Rings[i].MetCheck("Whilst infall computed");
		double r = Rings[i].Radius;
		double w = Rings[i].Width;
		origMass[i] = Rings[i].Gas.Mass();
		double sigma = PredictSurfaceDensity(r,w,newGas,Rd);
		double newMass = sigma * 2.0 * pi * r*w;
		perfectMasses[i] = newMass;
		perfectDeltas[i] = newMass - Rings[i].Gas.Mass();
		if (perfectDeltas[i] < 0)
		{
			perfect = false;
		}
		perf += newMass;
	}
	std::vector<double> realDeltas(Rings.size(),0.0);
	if (perfect)
	{		
		realDeltas = perfectDeltas;
	}
	else
	{
		realDeltas = IterativeFit(perfectDeltas,predictedInfall);
	}
	for (int i = 0; i < Rings.size(); ++i)
	{
		double target = origMass[i] + realDeltas[i];
		double required = target - Rings[i].Gas.Mass(); // reocmpute mass to account for mass dragged through disc
		InsertInfallingGas(i,required);	
		
		Rings[i].MetCheck("After Infall applied " + std::to_string(required));
		RingMasses[i] = Rings[i].Mass();
	}	
}

double mass_integrand(double x)
{
	return -exp(-x) * (x + 1);
}

double Galaxy::PredictSurfaceDensity(double radius, double width, double totalGasMass, double scaleLength)
{
	double r = radius;
	double w = width;
	double truncatedFactor =  (1 - (1 + Param.Galaxy.Radius/scaleLength)*exp(-Param.Galaxy.Radius/scaleLength)); // from the finite size of the galaxy -- not an infinite exponential!
	double prefactor = totalGasMass/(2 * pi * r * w)  / truncatedFactor;
	double upRadius = (r+w/2)/scaleLength;
	double downRadius = (r - w/2)/scaleLength;
	return prefactor * (mass_integrand(upRadius) - mass_integrand(downRadius));
}


void Galaxy::CompoundScattering(int currentTime, int timeStart, int timeEnd)
{
	for (int t= timeStart; t < timeEnd; ++t)
	{
		Migrator[t].Compound(Migrator[currentTime]);
	}
}
void Galaxy::ComputeScattering(int t)
{
	if (Param.Migration.DispersionOrder > 0)
	{
		Migrator[t].Create(RingMasses);
		
		LaunchParallelOperation(t,t,Compounding);
	}
	//~ std::cout << "Scattering Matrix at " << t << std::endl;
	//~ for (int i = 0; i < Rings.size(); ++i)
	//~ {
		//~ for (int j = 0; j < Rings.size(); ++j)
		//~ {
			//~ std::cout << std::setw(15) << Migrator[t].Grid[i][j];
		//~ }
		//~ std::cout << std::endl;
	//~ }
}

void Galaxy::ScatterYields(int time, int ringstart, int ringend)
{
	double absorbFrac = 1.0 - Param.Stellar.EjectionFraction;
	
	for (int t = 0; t <= time; ++t)
	{
			
			
		//grab migrator
		const std::vector<std::vector<double>> & migrator = Migrator[t].Grid;
		
		for (int i = ringstart; i < ringend; ++i)
			{
		
			
			//self absorb fraction
			double selfAbsorb = migrator[i][i];
			Rings[i].Gas.Absorb(Rings[i].Stars.YieldsFrom(t),absorbFrac*selfAbsorb);
			bool dispersionContinues = true;
			int distance = 1;
			double truncationCheck = 0.0;
			while (dispersionContinues)
			{
				int upGrab = i + distance;
				int downGrab = i - distance;
				double grabFractionUp = 0;
				double grabFractionDown = 0;
				if (upGrab < Rings.size())
				{
					grabFractionUp = migrator[i][upGrab];
					Rings[i].Gas.Absorb(Rings[upGrab].Stars.YieldsFrom(t),absorbFrac*grabFractionUp);
				}
				if (downGrab >= 0)
				{
					grabFractionDown = migrator[i][downGrab];
					Rings[i].Gas.Absorb(Rings[downGrab].Stars.YieldsFrom(t),absorbFrac*grabFractionDown);
				}
				truncationCheck = std::max(grabFractionUp,grabFractionDown);
				++distance;
				if (distance >= Rings.size() || truncationCheck < Param.Migration.DispersionTruncation)
				{
					dispersionContinues = false;
				}
			}
			
			
			
			if (Param.Galaxy.IGMAbsorbing.Value)
			{
				IGM.Absorb(Rings[i].Stars.YieldsFrom(t),1.0 - absorbFrac); //this step might be broken with the parallelisation....
			}
		}
		
	}
	
	for (int i =0; i < Rings.size(); ++i)
	{
		Rings[i].UpdateMemory(time);
	}

		
}

void Galaxy::ScatterGas(int time)
{
	const std::vector<std::vector<double>> & migrator = Migrator[time].Grid;
	int n = Rings.size();
	std::vector<double> oldMasses(n,0.0);
	for (int i = 0; i < n; ++i)
	{
		oldMasses[i] = Rings[i].Gas.ColdMass();
	}
	for (int i = 0; i < n ; ++i)
	{
		
		int lower = std::max(0,i - Param.Migration.DispersionOrder-1);
		int upper = std::min(n,i + Param.Migration.DispersionOrder+2);
		for (int j = lower; j < upper; ++j)
		{
			if (j != i)
			{
				double availableMass = oldMasses[j];
				double scatteredMass = availableMass * migrator[i][j];
				Rings[i].Gas.TransferColdFrom(Rings[j].Gas,scatteredMass);
			}
		}
	}
	
}

void Galaxy::SaveState(double t)
{
	SaveState_Mass(t);
	SaveState_Events(t);
	SaveState_Enrichment(t);
	Data.Log("\tSaved state at " + std::to_string(t) + "\n",3);
}
void Galaxy::SaveState_Mass(double t)
{
	std::stringstream output;
	if (t == 0)
	{
		output << MassHeaders() << "\n";
	}
	int tt  = round(t/Param.Meta.TimeStep);
	for (int i = 0; i < Rings.size(); ++i)
	{
		double Ms = Rings[i].Stars.AliveMass();
		double Mc = Rings[i].Gas.ColdMass();
		double Mh = Rings[i].Gas.HotMass();
		MassReport Mrr = Rings[i].Stars.DeadMass();
		double Mwd = Mrr.WD/1e9;
		double Mns = Mrr.NS/1e9;
		double Mbh = Mrr.BH/1e9;
		double Mr = Mrr.Total/1e9;
		double Mt = Ms + Mc + Mh + Mr;
		double Migm = IGM.Mass();
		std::vector<double> vals = {Rings[i].Radius, Rings[i].Area,Mt,Ms,Mc,Mh,Mwd,Mns,Mbh,Migm};
		output << t;
		for (int j = 0; j < vals.size(); ++j)
		{
			output << ", " << vals[j];
		}
		output << "\n";
	}	
	JSL::writeStringToFile(Param.Output.GalaxyMassFile,output.str());
	
}
std::string Galaxy::MassHeaders()
{
	return "Time, Radius, SurfaceArea, TotalMass, StellarMass, ColdGasMass, HotGasMass, WDMass, NSMass, BHMass,IGMMass";
}

void Galaxy::SaveState_Enrichment(double t)
{
	
	int tt = round(t / Param.Meta.TimeStep);
	
	//only save to file at simiulation end!
	if (tt == Param.Meta.SimulationSteps-1)
	{
		Data.UrgentLog("\tSaving Chemical makeup:    ");
		int bars = 0;
		std::stringstream outputAbsoluteCold;
		std::stringstream outputLogarithmicCold;
		std::stringstream outputAbsoluteHot;
		std::stringstream outputLogarithmicHot;
		for (int time = 0; time < Param.Meta.SimulationSteps; ++time)
		{
			Data.ProgressBar(bars,time,Param.Meta.SimulationSteps);
			for (int i  = 0; i < Rings.size(); ++i)
			{
				Rings[i].SaveChemicalHistory(time,outputAbsoluteCold,outputLogarithmicCold,outputAbsoluteHot,outputLogarithmicHot);
			}
		}
		
		
		JSL::writeStringToFile(Param.Output.AbsoluteColdGasFile,outputAbsoluteCold.str());
		JSL::writeStringToFile(Param.Output.LogarithmicColdGasFile,outputLogarithmicCold.str());
		JSL::writeStringToFile(Param.Output.AbsoluteHotGasFile,outputAbsoluteHot.str());
		JSL::writeStringToFile(Param.Output.LogarithmicHotGasFile,outputLogarithmicHot.str());
	}
}

void Galaxy::SaveState_Events(double t)
{
	
	int tt = round(t / Param.Meta.TimeStep);
	
	//only save to file at simiulation end!
	if (tt == Param.Meta.SimulationSteps-1)
	{
		Data.UrgentLog("\tSaving Stellar Event Rate: ");
		int bars = 0;
		std::stringstream eventOutput;
		for (int time = 0; time < Param.Meta.SimulationSteps-1; ++time)
		{
			Data.ProgressBar(bars,time,Param.Meta.SimulationSteps);
			for (int i  = 0; i < Rings.size(); ++i)
			{
				Rings[i].Stars.SaveEventRate(time,eventOutput);
			}
		}
		
		
		JSL::writeStringToFile(Param.Output.EventRateFile,eventOutput.str());
	}
}
