#include "NuclearDisk.h"

//#include "JSL.h"

//#include <string>

/*//Nuclear Disk Constructor in Work
NuclearDisk::NuclearDisk(InitialisedData & data): Data(data), Param(data.Param), IGM(GasReservoir::Primordial(data.Param.Galaxy.IGM_Mass,data.Param))
{
	galaxyDir = 'Test';

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
	double w = Param.Galaxy.RingWidth[0];
	double p = Param.Migration.MarkovDispersionStrength * Param.Meta.TimeStep / (w * w);
	Data.UrgentLog("\tMigration Matrices initialised. Characteristic transition probability: " + std::to_string(p) + "\n");
}
*/

// Main Evolution Loop
void NuclearDisk::Evolve()
{
	Data.UrgentLog("\tEvolution will occur between 0Gyr and " + std::to_string(Param.Meta.SimulationDuration) + "Gyr, across " + std::to_string(Param.Meta.SimulationSteps) + " steps.\n");
	double t = 0;

	int fullBar = Param.Meta.ProgressHashes;
	int currentBars = 0;

	Data.UrgentLog("\tStarting Nuclear Disk evolution: ");

	getBarInflow();
	Data.UrgentLog("after getBarInflow");

	int finalStep = Param.Meta.SimulationSteps - 1; // intentionally offset by 1!

	for (int timestep = 0; timestep < finalStep; ++timestep)
	{
		// std::cout << "Time " << timestep << std::endl;

		updateBarInflowResevoir(timestep);

		Infall(t, timestep);

		// std::cout << "Computing scattering" << std::endl;
		ComputeScattering(timestep);

		// std::cout << "Computing rngs" << std::endl;
		LaunchParallelOperation(timestep, Rings.size(), RingStep);

		// std::cout << "Computing scattering pt 2" << std::endl;
		if (timestep < finalStep)
		{
			LaunchParallelOperation(timestep, Rings.size(), Scattering);
			ScatterGas(timestep);
		}

		//~ std::cout << "Computing savestate" << std::endl;

		Data.ProgressBar(currentBars, timestep, finalStep);
		SaveState(t);
		t += Param.Meta.TimeStep;
	}
}

/*rewrite Infall so more gas is funneled onto Nuclear Disk than by classical onfall*/

void NuclearDisk::Infall(double t, int timestep)
{
	double pi = 3.141592654;
	double Rd = GasScaleLength(t);
	double oldGas = ColdGasMass();
	double predictedInfall = InfallMass(timestep);

	double newGas = oldGas + predictedInfall;

	std::vector<double> origMass(Rings.size(), 0.0);
	std::vector<double> perfectMasses(Rings.size(), 0.0);
	std::vector<double> perfectDeltas(Rings.size(), 0.0);
	bool perfect = true;
	double perf = 0;
	for (int i = 0; i < Rings.size(); ++i)
	{
		Rings[i].MetCheck("Whilst infall computed");
		double r = Rings[i].Radius;
		double w = Rings[i].Width;
		origMass[i] = Rings[i].Gas.ColdMass();
		double sigma = PredictSurfaceDensity(r, w, newGas, Rd);
		double newMass = sigma * 2.0 * pi * r * w;
		perfectMasses[i] = newMass;
		perfectDeltas[i] = newMass - Rings[i].Gas.ColdMass();
		if (perfectDeltas[i] < 0)
		{
			perfect = false;
		}
		perf += newMass;
	}
	std::vector<double> realDeltas(Rings.size(), 0.0);
	if (perfect)
	{
		realDeltas = perfectDeltas;
	}
	else
	{
		realDeltas = IterativeFit(perfectDeltas, predictedInfall);
	}
	for (int i = 0; i < Rings.size(); ++i)
	{
		double target = origMass[i] + realDeltas[i];
		double required = target - Rings[i].Gas.ColdMass(); // reocmpute mass to account for mass dragged through disc
		InsertInfallingGas(i, required);

		//~ Rings[i].MetCheck("After Infall applied " + std::to_string(required));
		RingMasses[i] = Rings[i].Mass();
	}
}

// accretion in the begining -> what to do with ring 0?
double NuclearDisk::InfallMass(int timestep)
{
	// qq
	//  if (timestep <= 100){
	//  	return 0;
	//  }
	double accretedInflow = IGM.ColdMass() * (1.0 - Param.NuclearDisk.ColdGasTransportLoss) + IGM.HotMass() * (1.0 - Param.NuclearDisk.HotGasTransportLoss);

	double coldbarmass = 0;
	for (int i = 0; i < ProcessCount; ++i)
	{
		for (int e = 0; e < ElementCount; ++e)
		{
			coldbarmass += hotBarInflow[timestep][i][e];
		}
	}

	// std::cout<<timestep << ' ' << accretedInflow  << ' ' <<accretedInflow/timestepsPerRing[timestep]<< ' ' << IGM.ColdMass()<< ' ' << IGM.HotMass() << ' '<< coldbarmass 	<<std::endl;
	return accretedInflow / timestepsPerRing[timestep];
}

/*put in: hotGasLoss, coldGasLoss
 */
void NuclearDisk::updateBarInflowResevoir(int timestep)
{
	IGM.Wipe();

	for (int i = 0; i < Rings.size(); ++i)
	{
		// IGM.Absorb(Rings[i].IGMBuffer);
		Rings[i].IGMBuffer.Wipe();
	}

	for (int p = 0; p < ProcessCount; ++p)
	{
		// for (int i = 0; i < ElementCount; ++i)
		// {
		// 	coldBarInflow[timestep][p][i] *= 10;
		// }

		// GasReservoir dummy = GasReservoir();

		Gas coldGas = Gas(coldBarInflow[timestep][p]);
		Gas hotGas = Gas(hotBarInflow[timestep][p]);

		SourceProcess source = (SourceProcess)p;

		GasStream processGasStream = GasStream(source, hotGas, coldGas);

		IGM[source].Absorb(processGasStream);

		// qq why does gas resevoir have a Paramsobject that is not initialised?
	}

	double elemmass = 0;
	for (int e = 0; e < ElementCount; ++e)
	{
		elemmass += coldBarInflow[timestep][0][e];
	}
	// std::cout << "IGM "<<IGM[(SourceProcess)(0)].ColdMass()<< ' ' << IGM.ColdMass()<< ' ' << elemmass<<std::endl;
}

// Reads in a line from the output files and converts them to double vectors representing the gas streams for different processes
std::vector<std::vector<double>> NuclearDisk::readAndSliceInput(std::vector<std::string> stringVector)
{

	std::vector<std::vector<double>> processVectors;

	stringVector.erase(stringVector.begin(), stringVector.begin() + 3);

	std::vector<double> doubleVector(stringVector.size());
	std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](std::string const &val)
				   { return std::stod(val); });

	// throw away first stream giving total abundance and then separate into vectors according to source processes
	for (int p = 0; p < ProcessCount; ++p)
	{
		processVectors.push_back(std::vector<double>(doubleVector.begin() + (p + 1) * ElementCount, doubleVector.begin() + (p + 2) * ElementCount));
	}

	return processVectors;
}

// linear growth of the bar with a initialisiation time period of rapid growth and subsequent slow growth
std::vector<int> NuclearDisk::barGrowthFunction()
{
	std::vector<int> barLengthInRings(Param.Meta.SimulationSteps);

	for (int timestep = 0; timestep < Param.Meta.SimulationSteps; ++timestep)
	{
		double barLength;
		double dtime = timestep * Param.Meta.TimeStep;

		// Data.UrgentLog(std::to_string(timestep) + ' ' + std::to_string(dtime) + ' ' + std::to_string(barLengthInRings.size()) +  '\n');

		if (dtime < Param.NuclearDisk.BarGrowthStart)
		{
			barLength = 0;
		}
		else if (dtime < Param.NuclearDisk.BarGrowthStart + Param.NuclearDisk.BarInitialiseTime)
		{
			barLength = (dtime - Param.NuclearDisk.BarGrowthStart) * (Param.NuclearDisk.BarInitialLength / Param.NuclearDisk.BarInitialiseTime);
		}
		else
		{
			double growthSpeed = (Param.NuclearDisk.BarFinalLength - Param.NuclearDisk.BarInitialLength) / (Param.Meta.SimulationDuration - Param.NuclearDisk.BarGrowthStart - Param.NuclearDisk.BarInitialiseTime);
			barLength = Param.NuclearDisk.BarInitialLength + (dtime - Param.NuclearDisk.BarGrowthStart - Param.NuclearDisk.BarInitialiseTime) * growthSpeed;
		}

		double ringwidth = Param.NuclearDisk.GalaxyRadius / Param.NuclearDisk.GalaxyRingCount;
		int ringToEmpty = (int)(barLength / ringwidth);

		barLengthInRings[timestep] = ringToEmpty;

		// std::cout<< barLength <<  ' ' << ringwidth << ' ' << ringToEmpty<< std::endl;
	}
	return barLengthInRings;
}

void NuclearDisk::checkTimeResolution(std::string galaxyFileCold, std::string galaxyFileHot){
	double t0 = 0;
	int i = 0;
	forLineVectorIn(
		galaxyFileCold, ', ',
		if (i > 0) 
		{
			// std::cout << FILE_LINE_VECTOR[0]<<std::endl;
			double t1 = std::stod(FILE_LINE_VECTOR[0]);
			// std::cout << t1-t0<< ' ' <<Param.Meta.TimeStep<< ' ' <<(t1-t0) - Param.Meta.TimeStep<<std::endl;
			if (t1 - t0 > 0.0)
			{
				if (std::abs(t1 - t0 - Param.Meta.TimeStep) > 0.0000001)
				{
					std::cout << "\n\t The chosen time resolution in the nuclear disk does not fit the input time resolution." << std::endl;
					exit(5);
				}
			}
			t0 = t1;
		} 
		++i;
	)

	t0 = 0;
	i = 0;
	forLineVectorIn(
		galaxyFileHot, ', ',
		if (i > 0) 
		{
			// std::cout << FILE_LINE_VECTOR[0]<<std::endl;
			double t1 = std::stod(FILE_LINE_VECTOR[0]);
			// std::cout << t1-t0<< ' ' <<Param.Meta.TimeStep<< ' ' <<(t1-t0) - Param.Meta.TimeStep<<std::endl;
			if (t1 - t0 > 0.0)
			{
				if (std::abs(t1 - t0 - Param.Meta.TimeStep) > 0.0000001)
				{
					std::cout << "\n\t The chosen time resolution in the nuclear disk does not fit the input time resolution." << std::endl;
					exit(5);
				}
			}
			t0 = t1;
		} 
		++i;
	)


}

void NuclearDisk::getBarInflow()
{
	std::string galaxyFileCold = "Output/" + galaxyDir + "/Enrichment_Absolute_ColdGas.dat";
	std::string galaxyFileHot = "Output/" + galaxyDir + "/Enrichment_Absolute_HotGas.dat";

	Data.UrgentLog(galaxyFileCold + '\n');

	// std::cout <<Param.Galaxy.Radius << " " << Param.Galaxy.RingCount << " " << Param.Galaxy.RingWidth[5] << "\n";

	std::vector<int> barLengthInRings = barGrowthFunction();

	timestepsPerRing.resize(Param.Meta.SimulationSteps);
	int sameRingCount = 1;
	for (int r = 1; r < barLengthInRings.size(); ++r)
	{
		// currently very last timestep doesn't accrete
		if ((barLengthInRings[r] == barLengthInRings[r - 1]) && (r != barLengthInRings.size() - 1))
		{
			++sameRingCount;
		}

		else if (r == barLengthInRings.size() - 1)
		{
			sameRingCount++;
			for (int incr = r - sameRingCount + 1; incr <= r; ++incr)
			{
				timestepsPerRing[incr] = sameRingCount;
			}
		}
		else
		{
			// std::cout << r - sameRingCount << ' ' << r << '\n';
			for (int incr = r - sameRingCount; incr < r; ++incr)
			{
				timestepsPerRing[incr] = sameRingCount;
			}
			sameRingCount = 1;
		}
	}

	// checking whether time resolutions match
	checkTimeResolution(galaxyFileCold, galaxyFileHot);

	int timestep = 0;

	forLineVectorIn(
		galaxyFileCold, ', ',
		std::string comp = std::to_string(barLengthInRings[timestep]) + ',';
		// std::string comp2 = std::to_string(timestep) + ',';

		if (FILE_LINE_VECTOR[1] == comp) {
			// std::cout<< timestep<<' ' <<FILE_LINE_VECTOR[1]<<std::endl;
			// std::cout << FILE_LINE<<std::endl;
			coldBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			++timestep;
		});

	// std::cout<<coldBarInflow[350]

	timestep = 0;
	forLineVectorIn(
		galaxyFileHot, ', ',
		std::string comp = std::to_string(barLengthInRings[timestep]) + ',';
		if (FILE_LINE_VECTOR[1] == comp) {
			hotBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			++timestep;
		});

	for (int i = 0; i < Param.Meta.SimulationSteps; i += 50)
	{
		Data.UrgentLog(std::to_string(i) + ' ');
		for (int e = 0; e < ElementCount; ++e)
		{
			Data.UrgentLog(std::to_string(coldBarInflow[i][0][e]) + ' ');
			// std::cout << i << ' '<< coldBarInflow[i][0][e];
		}
		// Data.UrgentLog( " tab ");
		// for (int e = 0; e < ElementCount; ++e)
		// {
		// 	Data.UrgentLog(std::to_string(coldBarInflow[i][1][e])+ ' ');
		// }
		Data.UrgentLog(" \n ");
	}
	Data.UrgentLog("\n");
	for (int i = 0; i < Param.Meta.SimulationSteps; i += 50)
	{
		for (int e = 0; e < ElementCount; ++e)
		{
			Data.UrgentLog(std::to_string(hotBarInflow[i][0][e]) + ' ');
			// std::cout << i << ' '<< coldBarInflow[i][0][e];
		}
		// Data.UrgentLog( " tab ");
		// for (int e = 0; e < ElementCount; ++e)
		// {
		// 	Data.UrgentLog(std::to_string(coldBarInflow[i][1][e])+ ' ');
		// }
		Data.UrgentLog(" \n ");
	}
}

double NuclearDisk::GasScaleLength(double t)
{

	double tdelay = Param.NuclearDisk.BarGrowthStart;
	double t0 = Param.Galaxy.ScaleLengthDelay;
	double tg = Param.Galaxy.ScaleLengthTimeScale;
	double tf = Param.Galaxy.ScaleLengthFinalTime;
	double R0 = Param.Galaxy.MinScaleLength;
	double Rf = Param.Galaxy.MaxScaleLength;

	if (t < tdelay)
	{
		// std::cout << t << ' ' << R0 <<  ' ' << Param.Meta.TimeStep <<'\n';
		return R0;
	}

	double N = 1.0 / (atan((tf - tdelay - t0) / tg) - atan(-t0 / tg));

	// std::cout << t << ' ' <<  R0 + N*(Rf - R0) * (atan( (t -tdelay - t0)/tg) - atan(-t0/tg)) << '\n';

	return R0 + N * (Rf - R0) * (atan((t - tdelay - t0) / tg) - atan(-t0 / tg));
}
