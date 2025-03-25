#include "GlobalParameters.h"

GlobalParameters::GlobalParameters()
{
	// anything here?
}

void GlobalParameters::Initialise(int argc, char *argv[])
{
	for (int i = 0; i < ParamMembers.size(); ++i)
	{
		ParamMembers[i]->Configure(argc, argv);
	}

	SaveInputs();

	for (int i = 0; i < ParamMembers.size(); ++i)
	{
		ParamMembers[i]->Initialise(Resources.ResourceRoot);
	}
}

void GlobalParameters::SaveInputs()
{

	std::string configOut = Output.Root.Value + "/" + Output.Config.Value;
	JSL::initialiseFile(configOut);

	std::stringstream output;
	for (int i = 0; i < ParamMembers.size(); ++i)
	{

		ParamMembers[i]->StreamContentsTo(output);
	}
	// JSL::initialiseFile(configOut);
	JSL::writeStringToFile(configOut, output.str());


	//~ write out massgrid 
	
	std::string MassGridFile = Output.Root.Value + "/Massgrid.dat";
	JSL::initialiseFile(MassGridFile);
	std::stringstream output ;
	output << "min_value,median_value,max_value\n";
	
	for (int i = 0; i < Stellar.MassResolution; ++i)
	{
		double x = Stellar.MassGrid[i];
		double w = Stellar.MassDeltas[i];

		output << x- 0.5*w << ","  << x<< "," << x + 0.5*w << "\n";
	}
	
	JSL::writeStringToFile(MassGridFile, output.str());
}
