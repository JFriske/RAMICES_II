#include "StellarPopulation.h"
#include <stdexcept>
IsoMass::IsoMass()
{
	MassIndex = 0;
	Count = 0;
	BirthIndex = 0;
	DeathIndex = 1e10;
}
IsoMass::IsoMass(int n, int m, double z, int birth, int death)
{
	MassIndex = m;
	Count = n;
	Metallicity = z;
	BirthIndex = birth;
	DeathIndex = death;
}




StellarPopulation::StellarPopulation(InitialisedData & data): Param(data.Param), IMF(data.IMF), SLF(data.SLF), CCSNYield(data.CCSNYield)
{
	Distribution.resize(Param.Stellar.MassResolution);
	internal_MassCounter = 0;
	IsDepleted = false;
	IsLifetimeMonotonic = true;
}

IsoMass & StellarPopulation::Relic()
{
	return ImmortalStars;
}

const IsoMass & StellarPopulation::Relic() const
{
	return ImmortalStars;
}

const IsoMass & StellarPopulation::operator [](int i) const
{
	if (i < Distribution.size() && i >= 0)
	{
		return Distribution[i];
	}
	else
	{
		throw std::runtime_error("You just tried to access a member of a Stellar Population which does not exist!");
	}
}
IsoMass & StellarPopulation::operator [](int i)
{
	if (i < Distribution.size() && i >= 0)
	{
		return Distribution[i];
	}
	else
	{
		throw std::runtime_error("You just tried to access a member of a Stellar Population which does not exist!");
	}
}

void StellarPopulation::FormStars(double formingMass, int timeIndex,double formingMetallicity)
{
	double NStarsFormed = IMF.FormationCount(formingMass);
	double budget = 0;
	
	int prevIndex = timeIndex;
	for (int i = Param.Stellar.MassResolution -1; i >= 0; --i)
	{
		//~ std::cout << "Mass " << i << "  " << formingMetallicity << std::endl;
		double m = Param.Stellar.MassGrid[i];
		double nStars = NStarsFormed * IMF.Weighting(i);
		budget +=  nStars * m/1e9;
		
		int deathIndex = timeIndex + SLF(i,formingMetallicity);
		Distribution[i] = IsoMass(nStars,i,formingMetallicity, timeIndex,deathIndex);
		
		//check for monotonicity
		if (deathIndex < prevIndex)
		{
			IsLifetimeMonotonic = false;
		}
		prevIndex = deathIndex;
	}
	//the remaining mass gets turned into immortal stars, al of which are assumed to have the minimum mortal mass
	double mInf = Param.Stellar.ImmortalMass;
	double effectiveImmortalCount = std::max(0.0,formingMass - budget)*1e9 / mInf;
	ImmortalStars = IsoMass(effectiveImmortalCount,mInf,formingMetallicity,timeIndex,1e10);
	
	internal_MassCounter += formingMass;
	DepletionIndex = Param.Stellar.MassResolution -1;
	
	if (std::isnan(internal_MassCounter))
	{
		std::cout << "Encountered critical errror " << std::endl;
	}
}

double StellarPopulation::Mass()
{
	
	if (std::isnan(internal_MassCounter))
	{
		std::cout << "Encountered critical errror " << std::endl;
		exit(5);
	}
	
	return internal_MassCounter;
}
bool StellarPopulation::Active()
{
	return !IsDepleted;
}
void StellarPopulation::Death(int time, GasReservoir & temporalYieldGrid, RemnantPopulation & remnants, GasReservoir & birthGas)
{
	if (IsLifetimeMonotonic)
	{
		MonotonicDeathScan(time, temporalYieldGrid, remnants,birthGas);
	}
	else
	{
		FullDeathScan(time);
	}
}
void StellarPopulation::MonotonicDeathScan(int time, GasReservoir & temporalYieldGrid, RemnantPopulation & remnants, GasReservoir & birthGas)
{

	while ( (Distribution[DepletionIndex].DeathIndex <= time || Distribution[DepletionIndex].Count == 0) && DepletionIndex >= 0)
	{
		//recover population information
		int nStars = Distribution[DepletionIndex].Count;
		int massID = Distribution[DepletionIndex].MassIndex;
		double starMass = Param.Stellar.MassGrid[massID];
		int birthID = Distribution[DepletionIndex].BirthIndex;
		double z = Distribution[DepletionIndex].Metallicity;
		 
		 
		 if (nStars > 0)
		 {
			double stellarMassReleased = nStars * starMass;
			Distribution[DepletionIndex].Count = 0;
			double gasMassReclaimed = stellarMassReleased / 1e9;
			internal_MassCounter -= gasMassReclaimed;
	
			if (starMass > Param.Yield.CCSN_MassCut)
			{
				CCSNYield(temporalYieldGrid,remnants,nStars,massID,z,birthID,birthGas);
			}
		}
		
		--DepletionIndex;
	}
	
	
	bool remainingStarsOutliveSimulation = Distribution[DepletionIndex].DeathIndex > Param.Meta.SimulationSteps;
	if (DepletionIndex < 0 || remainingStarsOutliveSimulation)
	{
		IsDepleted = true;
	}
	
	
	//~ GasStream output(CCSN,creation,Param.Thermal.HotInjection_CCSN);

}
void StellarPopulation::FullDeathScan(int time)
{
	throw std::runtime_error("You are calling death on a non-monotonic population, but this functionality does not exist yet!");
}
