#pragma once
#include "JSL.h"
#include "List.h"
#include "EnumSets.h"
using JSL::Argument;


//!The MetaValues contains variables associated with the base-level information about the sumulation - the number of cores to access, the timesteps etc. 
class MetaValues : public ParamList
{
	
	public:
		//! Controls whether the funky ASCII welcome message is played at the beginning of the code
		Argument<bool> WelcomeMessage = Argument<bool>(false, "welcome");

		//! The location of the funky ASCII welcome messgae
		Argument<std::string> WelcomeFile = Argument<std::string>("welcome.dat","welcome-file");
		
		//!The location of the directory which the code looks for its expected resource file structure
		Argument<std::string> ResourceRoot  = Argument<std::string>("Resources/","resource");
		
		//!The maximum number of parallel threads which can be active at any given time
		Argument<int> ParallelThreads = Argument<int>(1,"thread");
		
		//!The top level timestep used in the main chemical loop
		Argument<double> TimeStep  = Argument<double>(0.1,"timestep");
		
		//!The total duration of the chemical simulation
		Argument<double> SimulationDuration = Argument<double>(14.0,"duration");
		
		//!The number of timesteps in the simulation, computed from #SimulationDuration and #TimeStep
		int SimulationSteps;
		
		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array.
		MetaValues()
		{
			argPointers = {&WelcomeMessage, &WelcomeFile,&ResourceRoot,&ParallelThreads,&TimeStep,&SimulationDuration};
		};
		
		//! An overload of a normally empty function. Computes the value of #SimulationSteps
		virtual void Initialise(std::string resourceRoot);
};


class OutputValues : public ParamList
{
	public:
		//!The name of the output directory into which the output will be saved
		Argument<std::string> Root =  Argument<std::string>("Output/","output");
			
		//!The directory containing mass-type information
		Argument<std::string> GalacticDirectory = Argument<std::string>("Galactic/","galaxy-data");
		
		//!The name of the file containing galactic-scale mass information
		Argument<std::string> GalaxyMassFile = Argument<std::string>("Mass.dat","galaxy-mass-file");
		
		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array.
		OutputValues()
		{
			argPointers = {&Root, &GalacticDirectory,&GalaxyMassFile};
		};
		
		//! An overload of a normally empty function. Goes through and creates the necessary directory structure 
		virtual void Initialise(std::string resourceRoot);
};

//! The elemental suboptions contains variables and data associated with the solar abundances (and where to locate them), and how to extract and extrapolate the yield data from files. 
class ElementValues : public ParamList
{
	
	
	public:
		
		//! Human readable names for the elements, in the order associated with the ElementIDs. These names are primarily elemental symbols, except Metals, which uses "Z"
		std::vector<std::string> ElementNames;

		//! Solar abundances (in mass units) of the elements, in the order associated with the ElementIDs
		std::vector<double> SolarAbundances;
		
		//! The file in which the solar abundances can be found as a csv
		Argument<std::string> SolarAbundanceFile = Argument<std::string>("ChemicalData/SolarAbundances.dat","solar-values-file"); 

		//! The column of the solar abundance files which contains the ElementName for cross matching
		Argument<int> SolarAbundanceFileNameColumn = Argument<int>(0,"solar-values-name-col");
		
		//! The column of the solar abundance file which contains the relevant solar abundance value to be saved to memory
		Argument<int> SolarAbundanceFileDataColumn = Argument<int>(3,"solar-values-data-col");
		
		//! Enums to identify the theoretical basis for different yield tables
		enum SourceIDs {Orfeo,Marigo,Limongi,Maeder,Mixed};
		
		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array
		ElementValues()
		{
				argPointers = {&SolarAbundanceFile, &SolarAbundanceFileDataColumn, &SolarAbundanceFileNameColumn};
		};
		
		//! An overload of a normally empty function. Loads in the values fo the solar abundance data file into the SolarAbundances vector.
		virtual void Initialise(std::string resourceRoot);
		
		//! A fairly dumbly-written function which sorts the elemental symbols in ElementNames into the order specified by the global id-enum.
		void GiveElementsNames();
		
};

//! The subset of values associated with stars + their remnants
class StellarValues : public ParamList
{
	
	public:
		//!Minimum stellar mass that IMF can generate
		Argument<double> MaxStellarMass = Argument<double>(100,"max-mass");
		
		//!Maxmimum stellar mass that IMF can generate
		Argument<double> MinStellarMass = Argument<double>(0.5,"min-mass");
		
		//!Minimum Z that the ILM(??) can consider
		Argument<double> MinZ = Argument<double>(1e-7,"min-z");
		
		//!Maximum Z that the ILM(??) can consider
		Argument<double> MaxZ = Argument<double>(0.052,"max-z");

		//!The fraction of supernovae ejecta which is thrown into the IGM
		Argument<double> EjectionFraction = Argument<double>(0.45,"eject");
		
		//!For every 1 solar mass of stars which form, this fraction of gas is heated into the hot phase
		Argument<double> FeedbackFactor = Argument<double>(0,"mass-load");
		
		//! The fraction of white dwarfs which go SNIa
		Argument<double> SNIaFraction = Argument<double>(0.05,"sn1a-frac");
		
		//! The fraction of neutron stars which go NSM
		Argument<double> NSMFraction = Argument<double>(0.0001,"nsm-frac");
		
		//! The normal Kennicutt-Schmidt power law index
		Argument<double> SchmidtMainPower = Argument<double>(1.4,"schmidt-main");
		
		//! The low-density Kennicutt-Schmidt power law index
		Argument<double> SchmidtLowPower = Argument<double>(4.0,"schmidt-low");
		
		//! The density cut for the low/high density switchover in Schmidt power law
		Argument<double> SchmidtDensityCut = Argument<double>(4e-3,"schmidt-cut");
		
		//! The Schmidt prefactor
		Argument<double> SchmidtPrefactor = Argument<double>(1e-8,"schmidt-factor");
		
		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array
		StellarValues()
		{
			argPointers = {&MaxStellarMass, &MinStellarMass, &MinZ, &MaxZ, &EjectionFraction, &SNIaFraction, &NSMFraction,&SchmidtMainPower, &SchmidtLowPower, &SchmidtDensityCut, &SchmidtPrefactor};
		}
};

//!Thermal suboptions contain variables which deal with the thermal subroutines - cooling timescales injection fractions etc. 
class ThermalValues : public ParamList
{
	
	
	public:
		//!Fraction of CCSN ejecta which is put into the hot phase
		Argument<double> HotInjection_CCSN = Argument<double>(0.7,"fh-ccsn");
		
		//!Fraction of NSM ejecta which is put into the hot phase
		Argument<double> HotInjection_NSM = Argument<double>(0.4,"fh-nsm");
		
		//!Fraction of SNIa ejecta which is put into the hot phase
		Argument<double> HotInjection_SNIa = Argument<double>(0.99,"fh-sn1a");
		
		//! The exponential timescale over which the hot gas cools into the cold gas
		Argument<double> GasCoolingTimeScale = Argument<double>(1.0,"cool");
	
		
		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array
		ThermalValues()
		{
			argPointers = {&HotInjection_CCSN, &HotInjection_NSM, &HotInjection_SNIa, &GasCoolingTimeScale};
		}
};


//!The galaxy suboptions contians variables associated with the galaxy as a whole, such as the maximum radius, and various mass/infall properties
class GalaxyValues : public ParamList
{
	
	
	public:
		//! The number of annuli into which the galaxy is split
		Argument<int> RingCount = Argument<int>(100,"rings");
		
		
		//! The cutoff radius of the galaxy
		Argument<double> Radius = Argument<double>(20.0,"radius");
		
		//! Initial in-situ mass of the galaxy (assumed to be 100% gas)
		Argument<double> PrimordialMass = Argument<double>(8,"M0");
		
		
		//! Fraction of primordial gas which is hot
		Argument<double> PrimordialHotFraction = Argument<double>(0,"primordial-hot");
		
		//!Initial Mass of the IGM Reservoir
		Argument<double> IGM_Mass = Argument<double>(1e6,"igm-mass");
		
		//! The initial exponential scale length of the galaxy
		Argument<double> MinScaleLength = Argument<double>(1.75,"scale-length-min");
	
		//! The exponential scale length that the galaxy achieves at ScaleLengthFinalTime
		Argument<double> MaxScaleLength = Argument<double>(3.75,"scale-length-max");
		
		//! The delay time before the scale length begins to grow
		Argument<double> ScaleLengthDelay = Argument<double>(1.0,"scale-length-delay");
		
		//! The speed with which the scale length grows
		Argument<double> ScaleLengthTimeScale = Argument<double>(2.0,"scale-length-time");
		
		//! The time at which the scale length stops growing at becomes fixed
		Argument<double> ScaleLengthFinalTime = Argument<double>(12.0,"scale-length-final");

		//! The mass of the first (fast) exponential infall
		Argument<double> InfallMass1 = Argument<double>(4.5,"M1");
		
		//! The mass of the second (slow) exponential infall
		Argument<double> InfallMass2 = Argument<double>(45,"M2");
		
		//! The exponential timescale for the first (fast) exponential infall
		Argument<double> InfallTime1 = Argument<double>(0.3,"b1");
		
		//! The exponential timescale for the second (slow) exponential infall
		Argument<double> InfallTime2 = Argument<double>(14.0,"b2");

		//! A parameter to do with the inflow weighting scheme (icky)
		Argument<double> InflowParameterA = Argument<double>(0.33,"inflow-a");
		
		//! A parameter to do with the inflow weighting scheme (icky)
		Argument<double> InflowParameterB = Argument<double>(0.53,"inflow-b");
		
		//!maximum fraction which can be removed by SFR + associated feedback 
		Argument<double> MaxSFRFraction = Argument<double>(0.98,"max-sfr");

		//!Boring constructor -- slots in the relevant arguments into the ParamList::argPointer array
		GalaxyValues()
		{
			argPointers = {&RingCount, &PrimordialMass, &PrimordialHotFraction, &IGM_Mass, &Radius, &MinScaleLength, &MaxScaleLength, &ScaleLengthDelay, &ScaleLengthTimeScale, &ScaleLengthFinalTime, &InfallMass1, &InfallMass2, &InfallTime1, &InfallTime2, &InflowParameterA, &InflowParameterB, &MaxSFRFraction};
		}
	
};
