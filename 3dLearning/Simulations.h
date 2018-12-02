#pragma once

class Simulations
{
public:

	virtual void generateRequiredDistribution() = 0;
	virtual void gatherNecessaryData();
	virtual void fillSimulatedResult(int numberOfSim, int *resultArray)=0;

	virtual void consoleOutResult(int *result, int size);
};

class ArbitrarySimulations : public Simulations
{
public:
	~ArbitrarySimulations();
protected:
	// User entered or external input data
	float minValue, maxValue;
	int raimannSampleCount;// May be changed to constant


						   // Below data for calculations
	float midPoint;
	float minValueOffset, maxValueOffset;

	// Below data are calculated intermediates
	double* cdf;
public:

	virtual void gatherNecessaryData() override;
	virtual void fillSimulatedResult(int numberOfSim, int *resultArray) override;
};

class FunctionalSimulation : public Simulations
{
protected:

	virtual double getRandomVariable(double cdf)=0;
	virtual double yFunction(double x) { return x; };

public:
	virtual void generateRequiredDistribution() override {}
	virtual void gatherNecessaryData() override;
	virtual void fillSimulatedResult(int numberOfSim, int *resultArray) override;
};

class NormalDistributedArbitrarySim :public ArbitrarySimulations
{
private:

	const float mean = 0;
	const float sd = 1;

	inline double getProbablity(float x);

public:

	virtual void generateRequiredDistribution() override;
};

class ExponentDistSimulation : public FunctionalSimulation
{
protected:

	float rateOfDistribution = 0,rateOfAbsorbtion,rateOfScatter;

public:
	virtual double getRandomVariable(double cdf) override;

	float getAbsorptionRate() { return rateOfAbsorbtion; }

	float getScatterRate() { return rateOfScatter; }

	float getTransmitFallOffRate() { return rateOfDistribution; }

	virtual void gatherNecessaryData() override;
};

class BeerLambertLightPropagationSim : public ExponentDistSimulation
{

public:

	virtual void gatherNecessaryData() override;
};

class PhysicalLightPropagationSim : public FunctionalSimulation
{
private:


	int imageSize;
	unsigned int* lightPropagatedData;
	unsigned int* lightReflectedData;
	// Density of pixel per unit length of slab
	int pixelDensity;

	// Slab data
	float slabHeight, slabSize;
	// Using beer lambert propagation theory to find the random light propagation distribution
	BeerLambertLightPropagationSim* randomLightPropagationDistance;
	//-1 to 1 ,higher the value lesser the deflection will be
	float deflectionMag;

	// Photon data
	Vector3 enteringPoint=Vector3(), enteringDirection=Vector3(0,0,1);
	// Init energy of simulating photon group
	float photonGrpEnergy;

	int numberOfPasses;

	int nBounces, nAbsorbed, nTransmitted;

	float* reflectedDataAllPass;
	float* transmittedDataAllPass;

	int nCurrentPass;

	int rouletteSampleSize;

public:

private:

	void doSimulationPass(int nPhotons);

	void deflectPhoton(Vector3& direction, float g);

	float* smoothImageData(const unsigned int *imgData, int length, int width, int smoothSize);

protected:
	// Random Deflection of a photon group
	virtual double getRandomVariable(double g) override;

public:

	virtual ~PhysicalLightPropagationSim();

	virtual void gatherNecessaryData() override;

	virtual void fillSimulatedResult(int numberOfSim, int *resultArray) override;

	virtual void consoleOutResult(int *result, int size);
};