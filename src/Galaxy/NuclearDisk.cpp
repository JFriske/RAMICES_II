#include "NuclearDisk.h"
#include "Galaxy.h"
//double pi = 3.141592654;

//Main Evolution Loop
void NuclearDisk::Evolve()
{
	Data.UrgentLog("\tEvolution will occur between 0Gyr and " + std::to_string(Param.Meta.SimulationDuration) + "Gyr, across " + std::to_string(Param.Meta.SimulationSteps) + " steps.\n");
	double t = 0;


	int fullBar = Param.Meta.ProgressHashes;
	int currentBars = 0;
	
	Data.UrgentLog("\tStarting Nuclear Disk evolution: ");
	int finalStep = Param.Meta.SimulationSteps -1; // intentionally offset by 1!
	for (int timestep = 0; timestep < finalStep; ++timestep)
	{
		//~ std::cout << "Time " << timestep << std::endl
		getBarInflow();

		Infall(t);
		
		//~ std::cout << "Computing scattering" << std::endl;
		ComputeScattering(timestep);
		
		//~ std::cout << "Computing rngs" << std::endl;
		LaunchParallelOperation(timestep,Rings.size(),RingStep);
		
		//~ std::cout << "Computing scattering pt 2" << std::endl;
		if (timestep < finalStep)
		{
			LaunchParallelOperation(timestep,Rings.size(),Scattering);
			ScatterGas(timestep);
		}
		
		//~ std::cout << "Computing savestate" << std::endl;
		
		Data.ProgressBar(currentBars, timestep,finalStep);	
		SaveState(t);
		t += Param.Meta.TimeStep;
	}
}

void NuclearDisk::getBarInflow(){
    
}