#include "StarReservoir.h"

StarReservoir::StarReservoir(int parentRing, InitialisedData & data) : Data(data),Param(data.Param), ParentRing(parentRing), IMF(data.IMF), SLF(data.SLF), Remnants(data), YieldOutput(data.Param)
{
	StellarPopulation empty(Data,parentRing);
		
	for (int i = 0; i < Param.Meta.SimulationSteps+1; ++i)
	{
		Population.push_back(empty);
	}
	
	//Compute the parent ring surface area....needed for SFR so computing here is efficient!
	const double pi = 3.141592654;
	double width = Param.Galaxy.RingWidth[parentRing];
	double r = Param.Galaxy.RingRadius[parentRing];
	ParentArea = 2 * pi * r * width;
	Temp_Mass = 0;
	EventRate.resize(Param.Meta.SimulationSteps-1);
	PopulationIndex = 0;
}

double integratedSchmidt(double s0, double prefactor,double power, double t)
{
	if (power != 1)
	{
		double timeFactor = (power - 1)*prefactor* t;
		return pow( 1.0/pow(s0,power - 1) + timeFactor, 1.0/(1 - power));
	}
	else
	{
		return s0 * exp( - prefactor*t);
	}
}
double StarReservoir::SFR_GasLoss(double density)
{	
	double nBig = Param.Stellar.SchmidtMainPower;
	double nSmall = Param.Stellar.SchmidtLowPower;
	double sigmaCut = Param.Stellar.SchmidtDensityCut;
	double prefactor = Param.Stellar.SchmidtPrefactor * (1+Param.Stellar.FeedbackFactor);
	double power = nBig;
	if (density < sigmaCut)
	{
		power = nSmall;
		prefactor *= pow(sigmaCut,nBig - nSmall); // ensures the SFR is continuous
	}
	//integrates the SFR smoothly over the timestep (including Feedback losses), makes it impossible to losemore gas than you have
	
	double gasDensity = integratedSchmidt(density,prefactor,power,Param.Meta.TimeStep);
	
	return density - gasDensity;
	
}

void StarReservoir::Form(GasReservoir & gas)
{
	double z = gas.Metallicity();
	//~ std::cout << "I want to form stars" <<std::endl;
	double initMass = gas.ColdMass();
	//~ std::cout << "I have " << initMass << std::endl;
	//~ double initialTotalMass = gas.Mass() + Mass();
	
	////////   density version (old)
	double gasSurfaceDensity = gas.Mass() / ParentArea;
	double gasLossMass = std::max(0.0,ParentArea * SFR_GasLoss(gasSurfaceDensity));
	gasLossMass = std::min(gas.ColdMass() * 0.99, gasLossMass);
	//////// mass version (new)
	
	
	//~ double gasLossMass = SFR_GasLoss(gas.ColdMass());
	//~ std::cout << "I am about to lose " << gasLossMass << std::endl;
	double heatFrac = Param.Stellar.FeedbackFactor;
	
	
	double starMassFormed = 1.0/(1.0 + heatFrac) * gasLossMass;
	double feedbackMass = gasLossMass - starMassFormed;
	//~ std::cout << feedbackMass << "  " << starMassFormed << std::endl;
	gas.Deplete(starMassFormed,0.0);	
	//~ std::cout << "Attempting to heat the gas: " << gas.ColdMass() << "  " << gas.HotMass() << "  " << gas.Mass() << std::endl;
	gas.Heat(feedbackMass); 
	//~ std::cout << "After heating the gas: " << gas.ColdMass() << "  " << gas.HotMass() << "  " << gas.Mass() << std::endl;
	//~ std::cout << "Depleted successfully" <<std::endl;
	EventRate[PopulationIndex].StarMassFormed += starMassFormed;
	

	//~ std::cout<< "Calling form stars...." <<std::endl;
	int newStarCount = Population[PopulationIndex].FormStars(starMassFormed,PopulationIndex,z);
	//~ std::cout << newStarCount << "  " << PopulationIndex << std::endl;
	
	EventRate[PopulationIndex].NStarsFormed += newStarCount;
	//~ std::cout << "Checking gas mass" <<std::endl;
	if (gas.ColdMass() < 0)
	{
		std::cout << "ERROR!" << gas.ColdMass() << "  " << gas.HotMass() << std::endl;
		exit(3);
	}
	++PopulationIndex;
	//~ std::cout << "Done! " <<std::endl;
}

double StarReservoir::AliveMass()
{
	double m = 0;
	for (int i = 0; i < Population.size(); ++i)
	{
		m += Population[i].Mass();
	}
	return m;
}

void StarReservoir::PrintStatus(int t)
{
	std::string ringName = "Ring" + std::to_string(ParentRing) + "_stars.dat";
	std::stringstream output;
	int N = Param.Stellar.MassResolution;
	if (t == 0)
	{
		JSL::initialiseFile(ringName);
		output << "0, 0," << Param.Stellar.ImmortalMass-2;
		for (int j = 0; j < N; ++j)
		{
			output << ", " << Param.Stellar.MassGrid[j];
		}
		output << "\n";
	}
	
	
	for (int i = 0; i < t; ++i)
	{
		output << t << ", " << t-i << ", " <<Population[i].Relic().Count;
		for (int j = 0; j < N; ++j)
		{
			output << ", " << Population[i][j].Count;
		}
		output << "\n";
	}
	JSL::writeStringToFile(ringName, output.str());
}

void StarReservoir::Death(int currentTime, GasReservoir & birthGas)
{
	YieldOutput.WipeMemoryUpTo(currentTime);
	for (int i = 0; i < currentTime+1; ++i)
	{
		
		if (Population[i].Active())
		{
			Population[i].Death(currentTime, YieldOutput,Remnants, birthGas, EventRate[currentTime]);
		}
		
	}
	Remnants.Decay(currentTime,YieldOutput, EventRate[currentTime]);
}

const std::vector<GasStream> & StarReservoir::YieldsFrom(int t)
{
	//~ const std::vector<GasStream> output = YieldOutput.GetHistory(t);
	//~ std::cout << "Returning yield grid, XYZ incoming: " << std::endl;
	//~ for (int p = 0; p < ProcessCount; ++p)
	//~ {
		//~ std::cout << "Proc: " << output[p].Source << " X = " << output[p].Cold(Hydrogen) + output[p].Hot(Hydrogen) << "  Y = " << output[p].Cold(Helium) + output[p].Hot(Helium) << "  Z = " << output[p].Cold(Metals) + output[p].Hot(Metals) << " for total mass " << output[p].Mass() << std::endl;
	//~ }
	return YieldOutput.GetHistory(t);
}

MassReport StarReservoir::DeadMass()
{
	return Remnants.Mass();
}

void StarReservoir::SaveEventRate(int t, std::stringstream & output)
{
	if (t == 0 && ParentRing == 0)
	{
		EventRate[0].AddHeaders(output);
	}
	output << t * Param.Meta.TimeStep<< ", " << Param.Galaxy.RingRadius[ParentRing] << ", ";
	EventRate[t].Save(output,Param.Meta.TimeStep);
}

void StarReservoir::StealFrom(const StellarPopulation & mark, double fraction)
{
	StellarPopulation copy = mark;
	for (int i = 0; i < copy.Distribution.size(); ++i)
	{
		double n = copy.Distribution[i].Count * fraction;
		int intPart = n;
		
		double targetRoll = (n - intPart);
		double diceRoll = (double)rand() / RAND_MAX;
		if (diceRoll < targetRoll)
		{
			++intPart;
		}
		
		copy.Distribution[i].Count = intPart;
	}
	MigratedPopulation.push_back(copy);
}

void StarReservoir::Observations()
{
	std::string output = "";
	for (int i = 0; i < MigratedPopulation.size(); ++i)
	{
		std::vector<int> ms;
		
		double z = MigratedPopulation[i].Metallicity;
		double age = (Param.Meta.SimulationSteps - MigratedPopulation[i].BirthIndex) * Param.Meta.TimeStep; 
		
		for (int j = 0; j < MigratedPopulation[i].Distribution.size(); ++j)
		{
			if (MigratedPopulation[i].Distribution[j].Count > 0)
			{
				ms.push_back(MigratedPopulation[i].Distribution[j].MassIndex);
			}
		}
		std::vector<IsochroneEntry> data = Data.Isochrones.GetProperties(ms,z,age);
		std::vector<int> entryNumbers(ms.size(),0.0);
		for (int j = 0; j < ms.size(); ++j)
		{
			double entryFraction = 1e-3; //SOMETHING CLEVER HERE!
			
			double entryNumber = MigratedPopulation[i].Distribution[j].Count * entryFraction;
			int intPart = entryNumber;
		
			double targetRoll = (n - intPart);
			double diceRoll = (double)rand() / RAND_MAX;
			if (diceRoll < targetRoll)
			{
				++intPart;
			}
			
			
			
		}
		
		
	}
}
