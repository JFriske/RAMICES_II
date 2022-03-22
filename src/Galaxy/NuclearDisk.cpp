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

	Data.UrgentLog("before getBarInflow");
	getBarInflow();

	int finalStep = Param.Meta.SimulationSteps - 1; // intentionally offset by 1!

	for (int timestep = 0; timestep < finalStep; ++timestep)
	{
		//~ std::cout << "Time " << timestep << std::endl

		updateBarInflowResevoir(timestep);

		Infall(t);

		//~ std::cout << "Computing scattering" << std::endl;
		ComputeScattering(timestep);

		//~ std::cout << "Computing rngs" << std::endl;
		LaunchParallelOperation(timestep, Rings.size(), RingStep);

		//~ std::cout << "Computing scattering pt 2" << std::endl;
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

// /*rewrite Infall so more gas is funneled onto Nuclear Disk than by classical onfall*/
// void Infall(double t){
// 	std::cout <<"hi";
// }

/*take the first line from coldBarInflow and create a resevoir from that
take the first line from hot...
use those to do create the GasResevoir IGM
put in: hotGasLoss, coldGasLoss
delete the first lines
*/
void NuclearDisk::updateBarInflowResevoir(double timestep)
{

	for (int p = 0; p < ProcessCount; ++p)
	{	
		for(int i = 0; i< ElementCount;++i){
			coldBarInflow[timestep][p][i] *=10 ;
		}

		Gas coldGas = Gas(coldBarInflow[timestep][p]);
		Gas hotGas = Gas(hotBarInflow[timestep][p]);

		GasStream processGasStream = GasStream((SourceProcess)p, hotGas, coldGas);
		//GasReservoir barInflow = GasReservoir(param)
		IGM[(SourceProcess)p] = processGasStream;
		//qq why does gas resevoir have a Paramsobject that is not initialised?
	}
}

std::vector<std::vector<double>> NuclearDisk::readAndSliceInput(std::vector<std::string> stringVector)
{
	std::vector<std::vector<double>> processVectors;

	stringVector.erase(stringVector.begin(), stringVector.begin() + 3);

	std::vector<double> doubleVector(stringVector.size());
	std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](std::string const &val)
				   { return std::stod(val); });

	//throw away first stream giving total abundance and then separate into vectors according to source processes
	for (int p = 0; p < ProcessCount; ++p)
	{
		processVectors.push_back(std::vector<double>(doubleVector.begin() + (p+1) * ElementCount, doubleVector.begin() + (p + 2) * ElementCount));
	}

	return processVectors;
}

void NuclearDisk::getBarInflow()
{
	std::string galaxyFileCold = "Output/" + galaxyDir + "/Enrichment_Absolute_ColdGas.dat";
	std::string galaxyFileHot = "Output/" + galaxyDir + "/Enrichment_Absolute_HotGas.dat";

	Data.UrgentLog(galaxyFileCold+'\n');

	int ringToEmpty = 2;
	//insert function for ring to empty here


	forLineVectorIn(
		galaxyFileCold, ', ',
		std::string comp = std::to_string(ringToEmpty) + ',';
			if (FILE_LINE_VECTOR[1] == comp ) {
				coldBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			}
		);

	forLineVectorIn(
		galaxyFileHot, ', ',
		std::string comp = std::to_string(ringToEmpty) + ',';
			if (FILE_LINE_VECTOR[1] == comp ) {
				hotBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			}
		);


	for (int i = 0; i < 800; i +=100)
	{
		for (int e = 0; e < ElementCount; ++e)
		{
			Data.UrgentLog(std::to_string(coldBarInflow[i][0][e]) + ' ');
		}
		// Data.UrgentLog( " tab ");
		// for (int e = 0; e < ElementCount; ++e)
		// {
		// 	Data.UrgentLog(std::to_string(coldBarInflow[i][1][e])+ ' ');
		// }
		Data.UrgentLog( " \n ");
	}
}