#include "SimpleYield.h"


SimpleYield::SimpleYield(const GlobalParameters & param, YieldProcess yieldprocess): Param(param), Process(param.Yield.ProcessTypes[yieldprocess])
{
	if (param.Meta.Verbosity > 0)
	{
		std::cout << "\t" << Param.Yield.ProcessNames[Process] << " yield grid initialising...." << std::flush;
	}
	
	switch(yieldprocess)
	{
		case NSM:
		{
			NSM_Initialise();
			break;
		}
		case SNIa:
		{
			SNIa_Initialise();
			break;
		}
		default:
		{
			throw std::runtime_error("You have tried to initialise a yield grid for which there is no rule to create - ID = " +std::to_string(Process) + "...I am forced to quit");
		}
	}
	
	
	if (param.Meta.Verbosity > 0)
	{
		std::cout << "complete" << std::endl;
	}
}

void SimpleYield::operator()(GasReservoir & scatteringReservoir, int nObjects, int birthIndex) const
{
	RemnantInject(scatteringReservoir,nObjects,birthIndex);
}

void SimpleYield::RemnantInject(GasReservoir & scatteringReservoir, int nObjects, int birthIndex) const
{
	GasStream TempStream(Process);
	for (int e = 0; e < ElementCount; ++e)
	{
		ElementID elem = (ElementID)e;
		double amountInjected = nObjects * Grid[elem] / 1e9;
		TempStream.Cold(elem) = amountInjected * (1.0 - hotInjectionFraction);
		TempStream.Hot(elem) = amountInjected * hotInjectionFraction;
	}
	
	scatteringReservoir.AbsorbMemory(birthIndex,TempStream);
}
void SimpleYield::SNIa_Initialise()
{
	hotInjectionFraction = Param.Thermal.HotInjection_SNIa;
	
	Grid = std::vector<double>(ElementCount,0.0);
	
	Grid[Iron] = 0.8;
}
void SimpleYield::NSM_Initialise()
{
	hotInjectionFraction = Param.Thermal.HotInjection_NSM;
}
