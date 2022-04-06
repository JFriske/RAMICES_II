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
	 std::cout <<HotGasMass()<<' ';
	for (int r = 0; r < Param.Galaxy.RingCount; ++r)
	{
		double hotMass = Rings[r].Gas.HotMass();
		Rings[r].Gas.Deplete(0.0, Param.NuclearDisk.HotGasLossTimeStep * hotMass);
	}
	 std::cout<<HotGasMass()<<std::endl;
}

// accretion in the begining -> what to do with ring 0?
double NuclearDisk::InfallMass(int timestep)
{
	double accretedInflow = IGM.ColdMass() * (1.0 - Param.NuclearDisk.ColdGasTransportLoss) + IGM.HotMass() * (1.0 - Param.NuclearDisk.HotGasTransportLoss);

	double coldbarmass = 0;
	for (int i = 0; i < ProcessCount; ++i)
	{
		for (int e = 0; e < ElementCount; ++e)
		{
			coldbarmass += hotBarInflow[timestep][i][e];
		}
	}

	// std::cout<<timestep << ' ' << accretedInflow  << ' ' <<accretedInflow/timestepsPerRing[timestep]<< ' ' << IGM.ColdMass()<< ' ' << IGM.HotMass() << ' '<< coldbarmass 	<<std::endl;
	return accretedInflow / timestepsPerRing[timestep];
}

/*put in: hotGasLoss, coldGasLoss
 */
void NuclearDisk::updateBarInflowResevoir(int timestep)
{
	IGM.Wipe();

	// for (int i = 0; i < Rings.size(); ++i)
	// {
	// 	// IGM.Absorb(Rings[i].IGMBuffer);
	// 	Rings[i].IGMBuffer.Wipe();
	// }

	for (int p = 0; p < ProcessCount; ++p)
	{
		// std::cout << timestep << " process ";

		// std::cout << coldBarInflow[timestep][p][Magnesium] << ' ';

		Gas coldGas = Gas(coldBarInflow[timestep][p]);
		Gas hotGas = Gas(hotBarInflow[timestep][p]);

		SourceProcess source = (SourceProcess)p;

		GasStream processGasStream = GasStream(source, hotGas, coldGas);

		IGM[source].Absorb(processGasStream);
		// qq why does gas resevoir have a Paramsobject that is not initialised?
	}
	// std::cout << std::endl;

	double elemmass = 0;

	for (int e = 0; e < ElementCount; ++e)
	{
		elemmass += coldBarInflow[timestep][0][e];
	}

	//std::cout << timestep * 0.03 << " IGM " << IGM[(SourceProcess)(0)].ColdMass() << ' ' << IGM.ColdMass() << ' ' << elemmass << std::endl;
}

// Reads in a line from the output files and converts them to double vectors representing the gas streams for different processes
std::vector<std::vector<double>> NuclearDisk::readAndSliceInput(std::vector<std::string> stringVector)
{

	std::vector<std::vector<double>> processVectors;

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
std::vector<int> NuclearDisk::barGrowthFunction()
{
	std::vector<int> barLengthInRings(Param.Meta.SimulationSteps);

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
		int ringToEmpty = (int)(barLength / ringwidth);

		barLengthInRings[timestep] = ringToEmpty;

		// std::cout<< barLength <<  ' ' << ringwidth << ' ' << ringToEmpty<< std::endl;
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
	std::string galaxyFileCold = "Output/" + galaxyDir + "/Enrichment_Absolute_ColdGas.dat";
	std::string galaxyFileHot = "Output/" + galaxyDir + "/Enrichment_Absolute_HotGas.dat";

	Data.UrgentLog(galaxyFileCold + '\n');

	// std::cout <<Param.Galaxy.Radius << " " << Param.Galaxy.RingCount << " " << Param.Galaxy.RingWidth[5] << "\n";

	std::vector<int> barLengthInRings = barGrowthFunction();

	timestepsPerRing.resize(Param.Meta.SimulationSteps);
	int sameRingCount = 1;
	for (int r = 1; r < barLengthInRings.size(); ++r)
	{
		// currently very last timestep doesn't accrete
		if ((barLengthInRings[r] == barLengthInRings[r - 1]) && (r != barLengthInRings.size() - 1))
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
	}

	// checking whether time resolutions match
	checkTimeResolution(galaxyFileCold, galaxyFileHot);

	int timestep = 0;

	forLineVectorIn(
		galaxyFileCold, ', ',
		std::string comp = std::to_string(barLengthInRings[timestep]) + ',';
		if (FILE_LINE_VECTOR[1] == comp) {
			coldBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			++timestep;
		});

	timestep = 0;
	forLineVectorIn(
		galaxyFileHot, ', ',
		std::string comp = std::to_string(barLengthInRings[timestep]) + ',';
		if (FILE_LINE_VECTOR[1] == comp) {
			hotBarInflow.push_back(readAndSliceInput(FILE_LINE_VECTOR));
			++timestep;
		});

	// for (int i = 0; i < Param.Meta.SimulationSteps; i += 50)
	// {
	// 	Data.UrgentLog(std::to_string(i) + ' ');
	// 	for (int e = 0; e < ElementCount; ++e)
	// 	{
	// 		Data.UrgentLog(std::to_string(coldBarInflow[i][0][e]) + ' ');
	// 		// std::cout << i << ' '<< coldBarInflow[i][0][e];
	// 	}
	// 	Data.UrgentLog(" \n ");
	// }
	// Data.UrgentLog("\n");
	// for (int i = 0; i < Param.Meta.SimulationSteps; i += 50)
	// {
	// 	for (int e = 0; e < ElementCount; ++e)
	// 	{
	// 		Data.UrgentLog(std::to_string(hotBarInflow[i][0][e]) + ' ');
	// 		// std::cout << i << ' '<< coldBarInflow[i][0][e];
	// 	}
	// 	Data.UrgentLog(" \n ");
	// }
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
	double prefactor = 1 / (2 * pi * r * w);
	double upRadius = (r + w / 2) / scaleLength;
	double downRadius = (r - w / 2) / scaleLength;

	double total = (mass_integrand(upRadius) - mass_integrand(downRadius));

	double ringstrength = 0;

	double nuclearRingEdge = 2.0 * scaleLength + 0.5 * nuclearRingWidth;
	if (r > nuclearRingEdge)
	{
		total *= exp(-(r - nuclearRingEdge) / dropOffDelta);
	}

	total *= prefactor;

	if (r > 2.0 * scaleLength - 0.5 * nuclearRingWidth && r < nuclearRingEdge)
	{
		ringstrength = nuclearRingMassFraction / (1.0 - nuclearRingMassFraction) * (expNorm / (4.0 * pi * scaleLength * nuclearRingWidth));
		total += ringstrength;
	}

	double normFactor = expNorm + ringstrength * 4.0 * pi * scaleLength * nuclearRingWidth;

	return totalGasMass / normFactor * total;
}
