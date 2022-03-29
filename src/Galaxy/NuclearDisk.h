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
        std::string galaxyDir= "Test";


        void getBarInflow(); 
        void updateBarInflowResevoir(int timestep);
        double GasScaleLength(double t);
        void Infall(double time, int timestep);
        double InfallMass(int timestep);

        void checkTimeResolution(std::string galaxyFileCold, std::string  galaxyFileHot );
        std::vector<int> barGrowthFunction();


        std::vector<std::vector<double>> readAndSliceInput(std::vector<std::string> stringVector);

        std::vector<std::vector<std::vector<double>>> coldBarInflow;
        std::vector<std::vector<std::vector<double>>> hotBarInflow;

        std::vector<int> timestepsPerRing;
        
};
