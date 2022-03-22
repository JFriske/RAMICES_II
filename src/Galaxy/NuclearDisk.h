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
        void updateBarInflowResevoir(double timestep);
        //void Infall(double t);

        std::vector<std::vector<double>> readAndSliceInput(std::vector<std::string> stringVector);

        std::vector<std::vector<std::vector<double>>> coldBarInflow;
        std::vector<std::vector<std::vector<double>>> hotBarInflow;
        
};
