#pragma once
#include <vector>
#include "../Parameters/InitialisedData.h"
#include "IMF.h"
#include "../Gas/GasReservoir.h"
#include "RemnantPopulation.h"
#include "SLF.h"
#include "StarEvents.h"
//! A simple struct for tracking the number of stars of a given mass
class IsoMass
{
	public:
		int MassIndex;
		double Count;
		double Metallicity;
		int BirthIndex;
		int DeathIndex;
		IsoMass();
		IsoMass(int n, int m, double z, int birth, int death);
};


class StellarPopulation
{
	public:
		StellarPopulation(InitialisedData & data);
	
		void PrepareIMF();
		
		//!Returns the number of stars formed (spread across all mass grids)
		int FormStars(double formingMass, int timeIndex, double formingMetallicity);
		double Mass();
		IsoMass & Relic();
		const IsoMass & Relic() const;
		IsoMass & operator[](int i);
		const IsoMass & operator[](int i) const;
		bool Active();
		void Death(int time, GasReservoir & TemporalYieldGrid, RemnantPopulation & remnants, GasReservoir & birthGas, StarEvents & EventRate);
	private:
		const GlobalParameters & Param;
		IsoMass ImmortalStars;
		std::vector<IsoMass> Distribution;

		const IMF_Functor & IMF; 
		SLF_Functor & SLF;
		const YieldGrid & CCSNYield;
		
		bool IsLifetimeMonotonic;
		bool IsDepleted;
		int DepletionIndex;
		
		double internal_MassCounter;
		
		void MonotonicDeathScan(int time,GasReservoir & temporalYieldGrid, RemnantPopulation & remnants, GasReservoir & birthGas, StarEvents & eventRate);
		void FullDeathScan(int time);
		
		void RecoverMatter(int time,int nstars, int mass, GasReservoir & temporalYieldGrid, RemnantPopulation & remnants);
		
		Gas TempGas;
};
