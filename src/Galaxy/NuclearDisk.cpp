#include "NuclearDisk.h"

// Main Evolution Loop
void NuclearDisk::Evolve()
{
	Data.UrgentLog("\tEvolution will occur between 0Gyr and " + std::to_string(Param.Meta.SimulationDuration) + "Gyr, across " + std::to_string(Param.Meta.SimulationSteps) + " steps.\n");
	double t = 0;

	int fullBar = Param.Meta.ProgressHashes;
	int currentBars = 0;

	Data.UrgentLog("\tStarting Nuclear Disk evolution: ");

	getBarInflow();
	int finalStep = Param.Meta.SimulationSteps - 1; // intentionally offset by 1!

	for (int timestep = 0; timestep < finalStep; ++timestep)
	{
		// std::cout << "Time " << timestep << std::endl;

		updateBarInflowResevoir(timestep);
		SaveState_CGM(t,true);


		//		if (t < 9 || t > 9.5){
		//std::cout<<'here;'<<std::endl;
		Infall(t, timestep);
			//	}
		

		// std::cout << "Computing scattering" << std::endl;
		ComputeScattering(timestep);

		//std::cout<<"hatgasmass 1 "<<HotGasMass()<<std::endl;

		//  std::cout << "Computing rings" << std::endl;
		LaunchParallelOperation(timestep, Rings.size(), RingStep);

		//std::cout<<"hatgasmass 2 "<<HotGasMass()<<std::endl;
		//  std::cout << "Computing scattering pt 2" << std::endl;
		if (timestep < finalStep)
		{
			LaunchParallelOperation(timestep, Rings.size(), Scattering);
			ScatterGas(timestep);
		}

		LoseHotGas();

		//std::cout<<"hatgasmass end "<<HotGasMass()<<std::endl;

		//  std::cout << "Computing savestate" << std::endl;

		Data.ProgressBar(currentBars, timestep, finalStep);
		SaveState(t);

		SaveState_CGM(t,false);
		t += Param.Meta.TimeStep;
	}
}

/*rewrite Infall so more gas is funneled onto Nuclear Disk than by classical onfall*/

void NuclearDisk::Infall(double t, int timestep)
{
	double pi = 3.141592654;
	double Rd = GasScaleLength(t);
	double oldGas = ColdGasMass();
	double predictedInfall = InfallMass(timestep);

	// if (t >7.0){
	// 	predictedInfall = 0.0;
	// }

	// std::cout << t << ' ' << timestep << ' ' << predictedInfall << ' ' << oldGas << ' ' << predictedInfall / oldGas << std::endl;

	double newGas = oldGas + predictedInfall;

	std::vector<double> origMass(Rings.size(), 0.0);
	std::vector<double> perfectMasses(Rings.size(), 0.0);
	std::vector<double> perfectDeltas(Rings.size(), 0.0);
	bool perfect = true;
	double perf = 0;

	double exponentialDiskNorm = NormaliseSurfaceDensity(Rd);
	for (int i = 0; i < Rings.size(); ++i)
	{
		Rings[i].MetCheck("Whilst infall computed");
		double r = Rings[i].Radius;
		double w = Rings[i].Width;
		origMass[i] = Rings[i].Gas.ColdMass();
		double sigma = PredictSurfaceDensity(r, w, newGas, Rd, exponentialDiskNorm);
		double newMass = sigma * 2.0 * pi * r * w;
		perfectMasses[i] = newMass;
		perfectDeltas[i] = newMass - Rings[i].Gas.ColdMass();
		if (perfectDeltas[i] < 0)
		{
			perfect = false;
		}
		perf += newMass;
	}
	std::vector<double> realDeltas(Rings.size(), 0.0);
	if (perfect)
	{
		realDeltas = perfectDeltas;
	}
	else
	{
		realDeltas = IterativeFit(perfectDeltas, predictedInfall);
	}
	for (int i = 0; i < Rings.size(); ++i)
	{
		double target = origMass[i] + realDeltas[i];
		double required = target - Rings[i].Gas.ColdMass(); // reocmpute mass to account for mass dragged through disc
		InsertInfallingGas(i, required);

		//~ Rings[i].MetCheck("After Infall applied " + std::to_string(required));
		RingMasses[i] = Rings[i].Mass();
	}
}

void NuclearDisk::LoseHotGas()
{
	//std::cout <<HotGasMass()<<' ';
	for (int r = 0; r < Param.Galaxy.RingCount; ++r)
	{
		double hotMass = Rings[r].Gas.HotMass();
		Rings[r].Gas.Deplete(0.0, Param.NuclearDisk.HotGasLossTimeStep * hotMass);
	}
	 //std::cout<<HotGasMass()<<std::endl;
}

// accretion in the begining -> what to do with ring 0?
double NuclearDisk::InfallMass(int timestep)
{
	double accretedInflow = CGM.ColdMass() * (1.0 - Param.NuclearDisk.ColdGasTransportLoss) + CGM.HotMass() * (1.0 - Param.NuclearDisk.HotGasTransportLoss);

	// double coldbarmass = 0;
	// for (int i = 0; i < ProcessCount; ++i)
	// {
	// 	for (int e = 0; e < ElementCount; ++e)
	// 	{
	// 		coldbarmass += hotBarInflow[timestep][i][e];
	// 	}
	// }

	// std::cout<<timestep << ' ' << accretedInflow  << ' ' <<accretedInflow/timestepsPerRing[timestep]<< ' ' << CGM.ColdMass()<< ' ' << CGM.HotMass() << ' '<< coldbarmass 	<<std::endl;
	return accretedInflow / timestepsPerRing[timestep];
}



void NuclearDisk::SaveState_CGM(double t, bool early){
	std::stringstream outputAbsoluteCold_CGM;
	std::stringstream outputLogarithmicCold_CGM;
	std::stringstream outputAbsoluteHot_CGM;
	std::stringstream outputLogarithmicHot_CGM;

	int tt = round(t / Param.Meta.TimeStep);

	CGM_SaveChemicalHistory(tt, outputAbsoluteCold_CGM, outputLogarithmicCold_CGM, outputAbsoluteHot_CGM, outputLogarithmicHot_CGM);


	if(early)
	{
		JSL::writeStringToFile(Param.Output.LogarithmicCGMColdGasFileEarly, outputLogarithmicCold_CGM.str());
		JSL::writeStringToFile(Param.Output.LogarithmicCGMHotGasFileEarly, outputLogarithmicHot_CGM.str());
	}
	else
	{
		JSL::writeStringToFile(Param.Output.LogarithmicCGMColdGasFile, outputLogarithmicCold_CGM.str());
		JSL::writeStringToFile(Param.Output.LogarithmicCGMHotGasFile, outputLogarithmicHot_CGM.str());
	}
}

void NuclearDisk::CGM_SaveChemicalHistory(int t, std::stringstream & absoluteStreamCold, std::stringstream & logarithmicStreamCold, std::stringstream & absoluteStreamHot, std::stringstream & logarithmicStreamHot)
{	
	// std::cout<<"time " << t << " in CGM Log\n";

	std::vector<std::vector<double>> HotBuffer(ProcessCount + 1, std::vector<double>(ElementCount,0.0));
	std::vector<std::vector<double>> ColdBuffer(ProcessCount + 1, std::vector<double>(ElementCount,0.0));
		
	std::string basic = "";

	if (t == 0)
	{

		std::string headers = "Time";
		for (int p = -1; p < ProcessCount; ++p)
		{
			std::string processName;
			if (p > -1)
			{
					processName = Param.Yield.ProcessNames[p];
			}
			else
			{
				processName = "Total";
			}
			for (int e = 0; e < ElementCount; ++e)
			{
				std::string elementName = Param.Element.ElementNames[e];
			
			
				headers += ", " + processName+ "_" + elementName;
			}
		}
		basic = headers + ", ColdGasMass, HotGasMass, TotalMass, StepsPerRing \n";
		
	}


	basic += std::to_string(t*Param.Meta.TimeStep) ;
	
	absoluteStreamCold << basic;
	logarithmicStreamCold  << basic;
	absoluteStreamHot   << basic;
	logarithmicStreamHot   << basic;
	

	// std::cout<<"time " << t << "bef in CGM Log\n";

	const std::vector<GasStream> & target = CGM.Composition();

	// std::cout<<"time " << t << "after in CGM Log\n";
	
	double coldMass = 0;
	double hotMass = 0;
	for (int p = 0; p < ProcessCount; ++p)
	{
		double processCold = target[p].ColdMass();
		double processHot = target[p].HotMass();

		// std::cout<< processCold << " " << processHot<< "\n";
		
		coldMass += processCold;
		hotMass += processHot;
		for (int e = 0; e < ElementCount; ++e)
		{
			ElementID elem = (ElementID)e;
			double cold = target[p].Cold(elem);
			double hot = target[p].Hot(elem);
			
		// std::cout<< cold << " " << hot << " " << p << " "<< e<< "\n";

		// std::cout<< ColdBuffer[p][e]<<" "<< HotBuffer[p][e]<<"\n";

			if (p == 0)
			{
				// std::cout<< ColdBuffer[p][e]<<" "<< HotBuffer[p][e]<<"\n";
				ColdBuffer[p][e] = 0;
				HotBuffer[p][e] = 0;
			}
			// std::cout<< " here \n";
			ColdBuffer[0][e] += cold;
			HotBuffer[0][e] += hot;
			ColdBuffer[p+1][e] = cold/processCold;
			HotBuffer[p+1][e] = hot/processHot;
		}
	} 
	for (int e = 0; e < ElementCount; ++e)
	{
		ColdBuffer[0][e] /= (coldMass+1e-88);
		HotBuffer[0][e] /= (hotMass+1e-88);
	}
	

	for (int p = 0; p < ProcessCount + 1; ++p)
	{
		for (int e = 0; e < ElementCount; ++e)
		{
			double coldCorrect = coldMass;
			double hotCorrect = hotMass;
			if (p > 0)
			{
				coldCorrect = target[p-1].ColdMass();
				hotCorrect = target[p-1].HotMass();
			}


			neatLogAbs(ColdBuffer[p][e] * coldCorrect, absoluteStreamCold);
			neatLogAbs(HotBuffer[p][e] * hotCorrect, absoluteStreamHot);
					
			double logValueCold = log10(ColdBuffer[p][e] / Param.Element.SolarAbundances[e]);
			double logValueHot = log10(HotBuffer[p][e]/Param.Element.SolarAbundances[e]);
			
			neatLogLog(logValueCold,logarithmicStreamCold);
			neatLogLog(logValueHot,logarithmicStreamHot);
	
		}
	}


	absoluteStreamCold << ", "<< CGM.ColdMass() <<", " <<CGM.HotMass()<<", "<< CGM.Mass()<< ", " << timestepsPerRing[t] << "\n";
	logarithmicStreamCold  << ", "<< CGM.ColdMass()/ timestepsPerRing[t] <<", " <<CGM.HotMass()/ timestepsPerRing[t]<<", "<< CGM.Mass()/ timestepsPerRing[t]<< ", " << timestepsPerRing[t]<< "\n";
	// std::cout  << "rel "<< CGM.ColdMass()/ timestepsPerRing[t] <<", " <<CGM.HotMass()/ timestepsPerRing[t]<<", "<< CGM.Mass()/ timestepsPerRing[t]<< ", " << timestepsPerRing[t]<< "\n";
	// std::cout  << "abs "<< CGM.ColdMass()<<", " <<CGM.HotMass()<<", "<< CGM.Mass()<< ", " << timestepsPerRing[t]<< "\n";
	absoluteStreamHot   << ", "<< CGM.ColdMass() <<", " <<CGM.HotMass()<<", "<< CGM.Mass()<< ", " << timestepsPerRing[t] << "\n";
	logarithmicStreamHot  << ", "<< CGM.ColdMass() <<", " <<CGM.HotMass()<<", "<< CGM.Mass()<< ", " << timestepsPerRing[t]<< "\n";
}




/*put in: hotGasLoss, coldGasLoss
 */
void NuclearDisk::updateBarInflowResevoir(int timestep)
{
	CGM.Wipe();

	// for (int i = 0; i < Rings.size(); ++i)
	// {
	// 	// CGM.Absorb(Rings[i].CGMBuffer);
	// 	Rings[i].CGMBuffer.Wipe();
	// }

	for (int p = 0; p < ProcessCount; ++p)
	{
		// std::cout << timestep << " process ";

		//std::cout << coldBarInflow[timestep][p][Magnesium] << ' ';

		Gas coldGas = Gas(coldBarInflow[timestep][p]);
		Gas hotGas = Gas(hotBarInflow[timestep][p]);

		SourceProcess source = (SourceProcess)p;

		GasStream processGasStream = GasStream(source, hotGas, coldGas);

		CGM[source].Absorb(processGasStream);
		// qq why does gas resevoir have a Paramsobject that is not initialised?
	}
	// std::cout << std::endl;

	double elemmass = 0;

	for (int e = 0; e < ElementCount; ++e)
	{
		elemmass += coldBarInflow[timestep][0][e];
	}

	// std::cout << timestep * 0.01 << " CGM " << CGM[(SourceProcess)(0)].ColdMass() << ' ' << CGM.ColdMass() << ' ' << elemmass << std::endl;
}

// Reads in a line from the output files and converts them to double vectors representing the gas streams for different processes
std::vector<std::vector<double>> NuclearDisk::readAndSliceInput(std::vector<std::string> stringVector)
{

	std::vector<std::vector<double>> processVectors;

	// for (int i = 0; i < stringVector.size(); ++i){
	// 	std::cout<< stringVector[i]<<" ";
	// }
	// std::cout<<std::endl;

	stringVector.erase(stringVector.begin(), stringVector.begin() + 3);

	std::vector<double> doubleVector(stringVector.size());
	std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](std::string const &val)
				   { return std::stod(val); });

	// throw away first stream giving total abundance and then separate into vectors according to source processes
	for (int p = 0; p < ProcessCount; ++p)
	{
		processVectors.push_back(std::vector<double>(doubleVector.begin() + (p + 1) * ElementCount, doubleVector.begin() + (p + 2) * ElementCount));
	}
	return processVectors;
}

// linear growth of the bar with a initialisiation time period of rapid growth and subsequent slow growth
std::vector<double> NuclearDisk::barGrowthFunction()
{
	std::vector<double> barLengthInRings(Param.Meta.SimulationSteps);

	for (int timestep = 0; timestep < Param.Meta.SimulationSteps; ++timestep)
	{
		double barLength;
		double dtime = timestep * Param.Meta.TimeStep;

		// Data.UrgentLog(std::to_string(timestep) + ' ' + std::to_string(dtime) + ' ' + std::to_string(barLengthInRings.size()) +  '\n');

		if (dtime < Param.NuclearDisk.BarGrowthStart)
		{
			barLength = 0;
		}
		else if (dtime < Param.NuclearDisk.BarGrowthStart + Param.NuclearDisk.BarInitialiseTime)
		{
			barLength = (dtime - Param.NuclearDisk.BarGrowthStart) * (Param.NuclearDisk.BarInitialLength / Param.NuclearDisk.BarInitialiseTime);
		}
		else
		{
			double growthSpeed = (Param.NuclearDisk.BarFinalLength - Param.NuclearDisk.BarInitialLength) / (Param.Meta.SimulationDuration - Param.NuclearDisk.BarGrowthStart - Param.NuclearDisk.BarInitialiseTime);
			barLength = Param.NuclearDisk.BarInitialLength + (dtime - Param.NuclearDisk.BarGrowthStart - Param.NuclearDisk.BarInitialiseTime) * growthSpeed;
		}

		double ringwidth = Param.NuclearDisk.GalaxyRadius / Param.NuclearDisk.GalaxyRingCount;
		double ringToEmpty = (barLength / ringwidth);

		barLengthInRings[timestep] = ringToEmpty;

		 //std::cout<< barLength <<  ' ' << ringwidth << ' ' << ringToEmpty<< std::endl;
	}
	return barLengthInRings;
}

void NuclearDisk::checkTimeResolution(std::string galaxyFileCold, std::string galaxyFileHot)
{
	double t0 = 0;
	int i = 0;
	forLineVectorIn(
		galaxyFileCold, ', ',
		if (i > 0) {
			// std::cout << FILE_LINE_VECTOR[0]<<std::endl;
			double t1 = std::stod(FILE_LINE_VECTOR[0]);
			// std::cout << t1-t0<< ' ' <<Param.Meta.TimeStep<< ' ' <<(t1-t0) - Param.Meta.TimeStep<<std::endl;
			if (t1 - t0 > 0.0)
			{
				if (std::abs(t1 - t0 - Param.Meta.TimeStep) > 0.0000001)
				{
					std::cout << "\n\t The chosen time resolution in the nuclear disk does not fit the input time resolution." << std::endl;
					exit(5);
				}
			}
			t0 = t1;
		} ++i;)

		t0 = 0;
	i = 0;
	forLineVectorIn(
		galaxyFileHot, ', ',
		if (i > 0) {
			// std::cout << FILE_LINE_VECTOR[0]<<std::endl;
			double t1 = std::stod(FILE_LINE_VECTOR[0]);
			// std::cout << t1-t0<< ' ' <<Param.Meta.TimeStep<< ' ' <<(t1-t0) - Param.Meta.TimeStep<<std::endl;
			if (t1 - t0 > 0.0)
			{
				if (std::abs(t1 - t0 - Param.Meta.TimeStep) > 0.0000001)
				{
					std::cout << "\n\t The chosen time resolution in the nuclear disk does not fit the input time resolution." << std::endl;
					exit(5);
				}
			}
			t0 = t1;
		} ++i;)
}

void NuclearDisk::getBarInflow()
{
	std::string galaxyFileCold =  galaxyDir + "/Enrichment_Absolute_ColdGas.dat";
	std::string galaxyFileHot =   galaxyDir + "/Enrichment_Absolute_HotGas.dat";

	Data.UrgentLog(galaxyFileCold + '\n');

	// std::cout <<Param.Galaxy.Radius << " " << Param.Galaxy.RingCount << " " << Param.Galaxy.RingWidth[5] << "\n";

	std::vector<double> barLengthInRings = barGrowthFunction();

	timestepsPerRing.resize(Param.Meta.SimulationSteps);
	int sameRingCount = 1;
	for (int r = 1; r < barLengthInRings.size(); ++r)
	{
		//std::cout<< r << " " <<sameRingCount<<std::endl;
		// currently very last timestep doesn't accrete
		if (((int)barLengthInRings[r] == (int)barLengthInRings[r - 1]) && (r != barLengthInRings.size() - 1))
		{
			++sameRingCount;
			
		}

		else if (r == barLengthInRings.size() - 1)
		{
			sameRingCount++;
			for (int incr = r - sameRingCount + 1; incr <= r; ++incr)
			{
				timestepsPerRing[incr] = sameRingCount;
			}
		}
		else
		{
			// std::cout << r - sameRingCount << ' ' << r << '\n';
			for (int incr = r - sameRingCount; incr < r; ++incr)
			{
				timestepsPerRing[incr] = sameRingCount;
			}
			sameRingCount = 1;
		}

		//std::cout<< r << " end " <<sameRingCount<<std::endl;
	}

	//std::cout<<"HALLO/n";
    //Data.UrgentLog(galaxyFileCold + '\n');

	// checking whether time resolutions match
	checkTimeResolution(galaxyFileCold, galaxyFileHot);


	int timestep = 0;

	double barLength = barLengthInRings[timestep];

		std::vector<std::vector<double>> lowerRing;
		std::vector<std::vector<double>> upperRing;

	forLineVectorIn(
		galaxyFileCold, ', ',
		//Data.UrgentLog(galaxyFileCold + '\n');
		std::string compLower = std::to_string((int)barLength) + ',';
		std::string compUpper = std::to_string((int)(barLength)+1) + ',';
		if ( (int) barLength == Param.Galaxy.RingCount-1){
			compUpper = compLower;
		}

		if (FILE_LINE_VECTOR[1] == compLower) {
			lowerRing = readAndSliceInput(FILE_LINE_VECTOR);
			//std::cout<<lowerRing[1][1] <<'\n';
			//Data.UrgentLog(compLower +'\n');
		}
		if (FILE_LINE_VECTOR[1] == compUpper) {
			
			upperRing = readAndSliceInput(FILE_LINE_VECTOR);
			//Data.UrgentLog(compUpper +'\n');
			//Data.UrgentLog(upperRing[0][0] +'\n');
			std::vector<std::vector<double>> averageRing = lowerRing;
			for(int p = 0; p <= ProcessCount-1; ++p){
				//Data.UrgentLog(compUpper +'\n');
				for (int i = 0; i <= lowerRing[p].size(); ++i){
					double remainderFraction = fmod(barLength, (int)barLength);
					//std::cout <<"i "<< i<<" "<< p <<std::endl;
					if(barLength == 0){remainderFraction = 0;}
					//std::cout << barLength << " " << (int)barLength << " " <<remainderFraction<<" "<< (upperRing[p][0] - lowerRing[p][0] ) <<std::endl;
					//Data.UrgentLog(compUpper +'\n');
					averageRing[p][i] +=  remainderFraction *(upperRing[p][i] - lowerRing[p][i] );
				}
			}

			coldBarInflow.push_back(averageRing);
			++timestep;
		}
		);

	timestep = 0;
	forLineVectorIn(
		galaxyFileHot, ', ',
		//Data.UrgentLog(galaxyFileCold + '\n');
		std::string compLower = std::to_string((int)barLength) + ',';
		std::string compUpper = std::to_string((int)(barLength)+1) + ',';
		if ( (int) barLength == Param.Galaxy.RingCount-1){
			compUpper = compLower;
		}

		if (FILE_LINE_VECTOR[1] == compLower) {
			lowerRing = readAndSliceInput(FILE_LINE_VECTOR);
			//std::cout<<lowerRing[1][1] <<'\n';
			//Data.UrgentLog(compLower +'\n');
		}
		if (FILE_LINE_VECTOR[1] == compUpper) {
			
			upperRing = readAndSliceInput(FILE_LINE_VECTOR);
			//Data.UrgentLog(compUpper +'\n');
			//Data.UrgentLog(upperRing[0][0] +'\n');
			std::vector<std::vector<double>> averageRing = lowerRing;
			for(int p = 0; p <= ProcessCount-1; ++p){
				//Data.UrgentLog(compUpper +'\n');
				for (int i = 0; i <= lowerRing[p].size(); ++i){
					double remainderFraction = fmod(barLength, (int)barLength);
					//std::cout <<"i "<< i<<" "<< p <<std::endl;
					if(barLength == 0){remainderFraction = 0;}
					//std::cout << barLength << " " << (int)barLength << " " <<remainderFraction<<" "<< (upperRing[p][0] - lowerRing[p][0] ) <<std::endl;
					//Data.UrgentLog(compUpper +'\n');
					averageRing[p][i] +=  remainderFraction *(upperRing[p][i] - lowerRing[p][i] );
				}
			}

			hotBarInflow.push_back(averageRing);
			++timestep;
		}
		);

}

// change to more linear grow
double NuclearDisk::GasScaleLength(double t)
{

	double tdelay = Param.NuclearDisk.BarGrowthStart;
	double t0 = Param.Galaxy.ScaleLengthDelay;
	double tg = Param.Galaxy.ScaleLengthTimeScale;
	double tf = Param.Galaxy.ScaleLengthFinalTime;
	double R0 = Param.Galaxy.MinScaleLength;
	double Rf = Param.Galaxy.MaxScaleLength;
	double Ri = 0.0175 * Param.NuclearDisk.BarInitialLength;

	if (t < tdelay)
	{
		return R0;
	}
	else if (t < tdelay + t0)
	{
		return R0 + (t - tdelay) * (Ri - R0) / (t0);
	}
	else
	{
		return Ri + (t - tdelay - t0) * (Rf - Ri) / (Param.Meta.SimulationDuration - t0 - tdelay);
	}

	// double N = 1.0 / (atan((tf - tdelay - t0) / tg) - atan(-t0 / tg));

	// return R0 + N * (Rf - R0) * (atan((t - tdelay - t0) / tg) - atan(-t0 / tg));
}

double NuclearDisk::NormaliseSurfaceDensity(double scaleLength)
{
	double pi = 3.141592654;
	double nuclearRingWidth = 0.005; //change also below
	double dropOffDelta = 0.001;

	double sum = 0;

	for (int i = 0; i < Rings.size(); ++i)
	{
		double r = Rings[i].Radius;
		double w = Rings[i].Width;

		double upRadius = (r + w / 2) / scaleLength;
		double downRadius = (r - w / 2) / scaleLength;

		double total = (mass_integrand(upRadius) - mass_integrand(downRadius));

		double nuclearRingEdge = 2.0 * scaleLength + 0.5 * nuclearRingWidth;
		if (r > nuclearRingEdge)
		{
			total *= exp(-(r - nuclearRingEdge) / dropOffDelta);
		}
		sum += total;
	}
	return sum;
}

// adds a nuclear ring at 2*scaleLength with a mass fraction of the total mass and dropoff outside of it
double NuclearDisk::PredictSurfaceDensity(double radius, double width, double totalGasMass, double scaleLength, double expNorm)
{
	double pi = 3.141592654;
	double nuclearRingWidth = 0.005; //change also above
	double dropOffDelta = 0.001;
	double nuclearRingMassFraction = Param.NuclearDisk.NuclearRingMassFraction;
	double r = radius;
	double w = width;
	double prefactor = 1.0 / (2 * pi * r * w);
	double upRadius = (r + w / 2.0) / scaleLength;
	double downRadius = (r - w / 2.0) / scaleLength;

	double total = (mass_integrand(upRadius) - mass_integrand(downRadius));


	double nuclearRingEdge = 2.0 * scaleLength + 0.5 * nuclearRingWidth;
	if (r > nuclearRingEdge)
	{
		total *= exp(-(r - nuclearRingEdge) / dropOffDelta);
	}
	
	total *= 0.8 /expNorm;

	total *= prefactor;

	double stdev = 0.35*nuclearRingWidth;
	double mu = 2.0 * scaleLength;


	double nuclearRingGauss = 0.2/(r * 2.0 * pi * stdev * sqrt(2.0 * pi)) * exp( - 0.5 * (r-mu)*(r-mu)/(stdev*stdev));

	total += nuclearRingGauss;

	return totalGasMass * total;
}
