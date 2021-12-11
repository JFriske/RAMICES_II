#pragma once

#include <vector>
#include "../Parameters/GlobalParameters.h"
#include "../Gas/GasReservoir.h"
#include "StellarPopulation.h"
class StarReservoir
{
	public:
		StarReservoir(const GlobalParameters & Param, int parentRing);
		double Mass();
		
		void Birth(GasReservoir & gas);
		
		
	private:
		std::vector<StellarPopulation> Population;
		const GlobalParameters & Param;
		double SFR_GasLoss(double surfaceArea);
		const int ParentRing;
		double ParentArea;
		double Temp_Mass;
};
