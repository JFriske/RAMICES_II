#pragma once
#include "../Parameters/InitialisedData.h"
#include "Galaxy.h"

#include "Ring.h"
#include "../Gas/GasReservoir.h"

class NuclearDisk : public Galaxy
{
    public:
        //NuclearDisk(InitialisedData & Data);
        void Evolve();

    private:
        void getBarInflow(); 
        std::string galaxyFileName;
        
};
