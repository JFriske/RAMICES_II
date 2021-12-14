#include "ParameterLists.h"

void MetaValues::Initialise(std::string resourceRoot)
{
	SimulationSteps = ceil(SimulationDuration / TimeStep);
}

void OutputValues::Initialise(std::string resourceRoot)
{
	JSL::mkdir(Root.Value);
	GalacticDirectory.Value = Root.Value + "/" + GalacticDirectory.Value;
	JSL::mkdir(GalacticDirectory.Value);
	
	GalaxyMassFile.Value = GalacticDirectory.Value + "/" + GalaxyMassFile.Value;
	JSL::initialiseFile(GalaxyMassFile.Value);
	
	RingDirectory.Value = Root.Value + "/" + RingDirectory.Value;
	JSL::mkdir(RingDirectory.Value);
}


void ElementValues::GiveElementsNames()
{
	//pretty sure there must be a more intelligent way to do this....but I can guarantee that this preserves the ordering of the enum!
	ElementNames.resize(ElementCount);
	ElementNames[Hydrogen] = "H";
	ElementNames[Helium] = "He";
	ElementNames[Metals] = "Z";
	ElementNames[Iron] = "Fe";
	ElementNames[Oxygen] = "O";
	ElementNames[Magnesium] = "Mg";
	ElementNames[Carbon] = "C";
	ElementNames[Silicon] = "Si";
	ElementNames[Calcium] = "Ca";
	ElementNames[Manganese] = "Mn";
	ElementNames[Chromium]= "Cr";
	ElementNames[Cobalt] = "Co";
	ElementNames[Europium] = "Eu";
}

void ElementValues::Initialise(std::string resourceRoot)
{
	GiveElementsNames();
	
	//read in the solar abundance data from file
	std::string solarFile = resourceRoot + (std::string)SolarAbundanceFile;
	SolarAbundances = std::vector<double>(ElementCount,1e-20);
	int i = 0;

	forLineVectorIn(solarFile,',',
		if (i > 0)
		{
			int c = 0;
			while (c < ElementCount) 
			{
				if (FILE_LINE_VECTOR[SolarAbundanceFileNameColumn] == ElementNames[c])
				{
					SolarAbundances[c] = std::stod(FILE_LINE_VECTOR[SolarAbundanceFileDataColumn]);
					c = ElementCount;
				}
				++c;
			}
		}
		++i;
	);
}

double stepFraction(double targetMinStep, double width, int N)
{
	int nRaphson = 10;
	double x= 1;
	
	for (int i = 0; i < nRaphson; ++i)
	{
		double f;
		double fPrime;
		if (x == 1)
		{
			f= N - width;
			fPrime = N*(N-1)/2;
		}
		else
		{
			f= targetMinStep * (pow(x,N) - 1)/(x - 1) - width;
			fPrime = targetMinStep * ( (N-1)*pow(x,N) - N * pow(x,N-1) + 1)/pow(x -1,2);
		}
		x = x - f/fPrime;
	}
	return x;
}

void StellarValues::Initialise(std::string resourceRoot)
{
	//Allows for arbitrary mass steps without altering the rest of the code (hopefully!)
	MassGrid = std::vector<double>(MassResolution.Value);
	MassDeltas = std::vector<double>(MassResolution.Value);
	//~ double gridWidth = (MaxStellarMass - ImmortalMass)/MassResolution
	
	double minStepSize = 0.1;
	double alpha = stepFraction(minStepSize, MaxStellarMass - ImmortalMass,MassResolution.Value);
	std::cout << "Predicted alpha = " << alpha << std::endl;
	double sumFactor;
	if (alpha == 1)
	{
		sumFactor = MassResolution.Value;
	}
	else
	{
		sumFactor = (pow(alpha,MassResolution.Value) - 1.0)/(alpha - 1);
	}
	double w = (MaxStellarMass - ImmortalMass)/sumFactor;
	double x = ImmortalMass;
	
	for (int i = 0; i < MassResolution; ++i)
	{
		x += w/2;
		MassGrid[i] = x;
		MassDeltas[i] = w;
		x += w/2;
		w = w * alpha;
		std::cout << MassGrid[i] - MassDeltas[i]/2 << "->" << MassGrid[i] + MassDeltas[i]/2 << std::endl;
	}
	
	LogZGrid = std::vector<double>(LogZResolution);
	double zwidth = (MaxLogZ - MinLogZ)/(LogZResolution-1);
	for (int i =0; i < LogZResolution; ++i)
	{
		LogZGrid[i] = MinLogZ + (i)*zwidth;
	}
	LogZDelta = zwidth;
}
