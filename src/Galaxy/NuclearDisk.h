#pragma once
#include "../Parameters/InitialisedData.h"
#include "Galaxy.h"

#include "Ring.h"
#include "../Gas/GasReservoir.h"

#include "JSL.h"


class NuclearDisk : public Galaxy
{
    public:
        //NuclearDisk(InitialisedData & Data);
        using Galaxy::Galaxy;

        void Evolve();

    private:
        std::string galaxyDir= Param.NuclearDisk.GlobalReadInDir;


        void getBarInflow(); 
        void updateBarInflowResevoir(int timestep);
        double GasScaleLength(double t);
        double NormaliseSurfaceDensity(double scaleLength);
        double PredictSurfaceDensity(double radius, double width, double totalGasMass, double scaleLength, double expNorm);
        double PredictSurfaceDensityInflow(double radius, double width, double totalGasMass, double scaleLength, double expNorm, double time);
        double PredictSurfaceDensityOnfall(double radius, double width, double totalGasMass, double scaleLength, double expNorm, double time);
        void Infall(double time, int timestep);
        double InfallMass(int timestep);

        void InsertInfallingGasInflow(int ring, double amount, double scalelength);
        void InsertInfallingGasOnfall(int ring, double amount, double scalelength);


        void LoseHotGas();

        void CGM_SaveChemicalHistory(int t, std::stringstream & absoluteStreamCold, std::stringstream & logarithmicStreamCold, std::stringstream & absoluteStreamHot, std::stringstream & logarithmicStreamHot);
        
        void SaveState_CGM(double t, bool early);

        void checkTimeResolution(std::string galaxyFileCold, std::string  galaxyFileHot );
        std::vector<double> barGrowthFunction();


        std::vector<std::vector<double>> readAndSliceInput(std::vector<std::string> stringVector);

        std::vector<std::vector<std::vector<double>>> coldBarInflow;
        std::vector<std::vector<std::vector<double>>> hotBarInflow;

        std::vector<int> timestepsPerRing;
        
};
