#pragma once

//~ #include "../Stars/RemnantPopulation.h"
#include "../Parameters/GlobalParameters.h"
#include "../Gas/GasReservoir.h"

//~ #include <ostream>

struct RemnantOutput
{
	RemnantType Type;
	double Mass;
};

class YieldGrid
{
	public:
		const SourceProcess Process;
		YieldGrid(const GlobalParameters & param, SourceProcess Process);
		
		RemnantOutput operator()(GasReservoir & scatteringReservoir, int Nstars, int mass, double z, int birthIndex, GasReservoir & birthReservoir) const;
		
	private:
		const GlobalParameters & Param;
		std::vector<std::vector<std::vector<double>>> Grid;

		double hotInjectionFraction;
		
		void CCSN_Initialise();
		
		void AGB_Initialise();
		
		void InitialiseLargeGrid(int mSize, int zSize);
		//~ GasStream TempStream;
		//allows the grid size to be truncated for CCSN etc.
		int MassOffset;
		
		RemnantOutput StellarInject( GasReservoir & scatteringReservoir, int Nstars, int mass, double z, int birthIndex, GasReservoir & birthReservoir) const;
		
};
