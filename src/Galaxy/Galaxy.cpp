#include "Galaxy.h"

Galaxy::Galaxy(Options * opts)
{
	log(1) << "\nMain Galaxy object initialised\n";
	Opts = opts;
	
	double dr = Opts->Galaxy.MaxRadius / Opts->Simulation.NRings;
	
	//initialise galaxy masses

	Mass.ColdGas = 0.0;
	Mass.HotGas = 0.0;
	Mass.Stellar = 0.0;
	Mass.Remnant = 0.0;
	
	for (int i = 0; i < Opts->Simulation.NRings; ++i)
	{
		Ring R = Ring(opts,i,dr,this);
		Rings.push_back(R);
		
		Mass.ColdGas += R.Gas.ColdMass;
		Mass.HotGas += R.Gas.HotMass;
	}
	GasExcess = 0;
	IGM = GasReservoir(Opts);
	double igmMass = 100.0;
	IGM.SetPrimordial(igmMass);
	
	
	//Demosthenes = YieldGrid(opts);
	

}

double Galaxy::InfallRate(double t)
{
	double sum = 0.0;
	
	for (int i = 0; i < Opts->Galaxy.InfallMasses.size(); ++i)
	{
		double Mi = Opts->Galaxy.InfallMasses[i];
		double Bi = Opts->Galaxy.InfallTimeScales[i];
		
		sum += Mi/Bi * exp(-t/Bi);
	}

	return sum;
}

double Galaxy::GasScaleLength(double t)
{
	double t0 = Opts->Galaxy.ScaleLengthDelay;
	double tg = Opts->Galaxy.ScaleLengthGrowth;
	double tf = Opts->Galaxy.ScaleLengthTimeTether;
	
	double N = 1.0/(atan((tf - t0)/tg) - atan(-t0/tg));
	
	double R0 = Opts->Galaxy.MinScaleLength;
	double Rf = Opts->Galaxy.MaxScaleLength;
	
	return R0 + N*(Rf - R0) * (atan( (t - t0)/tg) - atan(-t0/tg));

}

void Galaxy::UpdateGasMass(double t)
{
	double dt = Opts->Simulation.TimeStep;
	int nRings = Opts->Simulation.NRings;
	double totalInfall = InfallRate(t)*dt;

	double newR = GasScaleLength(t);
	

	Mass.ColdGas = 0.0;
	Mass.HotGas = 0.0;
	for (int ringID = 0; ringID < nRings; ++ringID)
	{
		Mass.ColdGas += Rings[ringID].Gas.ColdMass;
		Mass.HotGas += Rings[ringID].Gas.HotMass;
	}
	double newMass = (Mass.ColdGas + Mass.HotGas + totalInfall - GasExcess);

	//actual newMass might differ from requested, so keep running sum, just in case
	
	double coldGas = 0;
	double hotGas = 0;
	for (int ringID = 0; ringID < nRings; ++ringID)
	{
		GasRequest request = Rings[ringID].AccretionRequest(t,newMass, newR);
		
		GasReservoir igm = PullIGM(request.IGM);
		
		igm.AddTo(&Rings[ringID].Gas);
		
		if (ringID < nRings - 1)
		{
			Rings[ringID].Gas.TakeFrom(&Rings[ringID+1].Gas,request.Disc);
		}


		Rings[ringID].UpdateInternalProperties();

		
		coldGas += Rings[ringID].Gas.ColdMass;
		hotGas  += Rings[ringID].Gas.HotMass;
	}
	
	Mass.ColdGas = coldGas;
	Mass.HotGas = hotGas;
	
	//GasExcess prevents the weird behaviour due to the inside-out protocols in the inner disk from causing newgas to be continually accreted
	GasExcess = std::max(0.0, Mass.ColdGas + Mass.HotGas - newMass);
}

GasReservoir Galaxy::PullIGM(double m)
{
	//simply copies the current IGM abundance into a new object and sets it to 100% cold gas. Does not alter the IGM.
	GasReservoir subset = IGM;
	subset.Mass = m;
	subset.ColdMass = m;
	subset.HotMass = 0;
	
	return subset;
}

void Galaxy::PushIGM(int donorID, double coldMass, double hotMass)
{
	Rings[donorID].Gas.GiveTo(&IGM, coldMass, hotMass);
}


template<class T>
void PushVectorLineToFile(std::fstream * file, std::vector<T> const & vec)
{
	int width = 12;
	for (T entry : vec)
	{
		*file << std::setw(width) << entry << ",";
	}
	*file << "\n";
}

void Galaxy::OpenLogs()
{
	std::string root = Opts->Simulation.FileRoot;
	
	RingState.open(root + Opts->Simulation.RingStateFile,std::fstream::out);
	std::vector<std::string> ringHeaders = Rings[0].PropertyHeaders();
	
	PushVectorLineToFile(&RingState, ringHeaders);
	
	GalaxyState.open(root + Opts->Simulation.GalaxyStateFile,std::fstream::out);
	std::vector<std::string> galaxyHeader = PropertyHeaders();
	
	PushVectorLineToFile(&GalaxyState, galaxyHeader);
}

void Galaxy::CloseLogs()
{
	RingState.close();
	GalaxyState.close();
}

void Galaxy::SaveState(double t)
{
	log(3) << "\tsaving simstate\n";

	std::vector<double> galaxyProperties = ReportProperties(t);
	PushVectorLineToFile(&GalaxyState, galaxyProperties);
	
	for (int id = 0; id < Opts->Simulation.NRings; ++id)
	{
		std::vector<double> ringProperties = Rings[id].ReportProperties(t);
		PushVectorLineToFile(&RingState, ringProperties);
	}
}

std::vector<std::string> Galaxy::PropertyHeaders()
{
	std::vector<std::string> headers = {"Time", "GasScaleLength", "TotalMass","ColdMass","HotMass","StellarMass","RemnantMass"};
	return headers;
}
std::vector<double> Galaxy::ReportProperties(double t)
{
	std::vector<double> properties = {t,GasScaleLength(t), Mass.Total(), Mass.ColdGas, Mass.HotGas, Mass.Stellar, Mass.Remnant};
	return properties;
}

void Galaxy::Evolve()
{
	OpenLogs();
	
	double dt = Opts->Simulation.TimeStep;
	int timeStep = 0;
	int nRings = Opts->Simulation.NRings;
	
	
	for (double t = 0; t < Opts->Simulation.FinalTime; t+=dt)
	{
		log(2) << "Simulation timestep t = " + std::to_string(t);
		SaveState(t);
		
		ChemicalEvolution(t);
		
		
		++timeStep;
	}
	
	SaveState(Opts->Simulation.FinalTime+dt);
	
	CloseLogs();
}


void Galaxy::ChemicalEvolution(double t)
{
	int nRings = Opts->Simulation.NRings;
	Mass.Stellar = 0;
	for (int i = 0; i < nRings; ++i)
	{
		Rings[i].FormStars(t);
		
		Mass.Stellar += Rings[i].Stars.Mass();
	}
	
	UpdateGasMass(t);
}


