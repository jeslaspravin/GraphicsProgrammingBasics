#include "stdafx.h"
#include "Simulations.h"
#include <cmath>
#include <algorithm>
#include <map>

void Simulations::consoleOutResult(int *result, int size)
{
	static const string fileOuput = "C:/Users/JesJas/Documents/Visual Studio 2017/CppProjects/3dLearning/PlotData/SimData.dt";

	ofstream fos(fileOuput);
	
	map<int, int> resultCount;
	for (int i = 0; i<size; i++)
	{
		resultCount[result[i]]++;
	}

	float multiplier = size / (float)1000;
	for (int i = 0; i < 50; i++) cout << "-";
	cout << '\n';
	for (auto pairItem = resultCount.begin(); pairItem != resultCount.end(); ++pairItem)
	{
		fos << pairItem->first << "  " << pairItem->second<<"\n";
		cout << pairItem->first << ((pairItem->second>0) ? "\t|[" : "\t|");
		for (float i = 0; i < pairItem->second; i += 5 * multiplier)
		{
			cout << "|";
		}
		cout << ((pairItem->second > 0) ? "]" : "") << '\n';
	}
	for (int i = 0; i < 50; i++) cout << "-";
	cout << "\n";
}


void Simulations::gatherNecessaryData()
{
#if WIN32 | _WIN32
	srand(1117);
#else
	srand48(1117);
#endif
}

ArbitrarySimulations::~ArbitrarySimulations()
{
	delete[] cdf;
}

void ArbitrarySimulations::gatherNecessaryData()
{
	Simulations::gatherNecessaryData();
	GET_STREAM_DATA("Enter number of discrete data to generate from sample range : ", raimannSampleCount);
	GET_TWO_STREAM_DATA("Enter Minimum and Maximum values of simulation : ", minValue, maxValue);

	cdf = new double[raimannSampleCount + 1];
	cdf[0] = 0;
	cdf[raimannSampleCount] = 1;

	if (minValue > maxValue)
	{
		cout << "Max and Min value out of order swapping them";
		maxValue = (minValue + maxValue) - (minValue = maxValue);// Swapping values
	}

	midPoint = (maxValue + minValue) / 2;
	minValueOffset = minValue - midPoint;
	maxValueOffset = maxValue - midPoint;
}

void ArbitrarySimulations::fillSimulatedResult(int numberOfSim, int *resultArray)
{
	for(int simCount=0;simCount<numberOfSim;simCount++)
	{
#if WIN32 | _WIN32
		float uniformRandom = float(rand()) / RAND_MAX;
#else
		float uniformRandom = drand48();
#endif

		double* selectedCDF = lower_bound(cdf, cdf + raimannSampleCount, uniformRandom);//Gives array element ptr which is <= value
		int minCDFIndex = selectedCDF-cdf- 1;
		int maxCDFIndex = minCDFIndex + 1;

		float t = float((uniformRandom - cdf[minCDFIndex]) / (cdf[maxCDFIndex] - cdf[minCDFIndex]));

		float offsetValue = (minCDFIndex + t) / raimannSampleCount * (maxValueOffset - minValueOffset);
		offsetValue += minValueOffset;

		if (resultArray + simCount)
		{
			resultArray[simCount] = (int)(midPoint + offsetValue);
		}
		else
		{
			throw "Array Index out of bound for result Array";
		}
	}	
}


inline double NormalDistributedArbitrarySim::getProbablity(float x)
{
	return exp(-1 * (x - mean)*(x - mean) / (2 * sd*sd)) / (sd * sqrt(2 * M_PI));
}

void NormalDistributedArbitrarySim::generateRequiredDistribution()
{
	float dx = (maxValueOffset - minValueOffset) / raimannSampleCount;
	for (int count = 1; count < raimannSampleCount; count++)
	{
		// Mid Raimann Sum Method
		float x = minValueOffset + (count + 0.5f) / (float)raimannSampleCount * (maxValueOffset - minValueOffset);
		double pr = getProbablity(x)*dx;
		cdf[count] = cdf[count - 1] + pr;
		printf("CDF of %0.2f is %0.3f\n", x, cdf[count]);
	}
}

void FunctionalSimulation::gatherNecessaryData()
{
	Simulations::gatherNecessaryData();
}

void FunctionalSimulation::fillSimulatedResult(int numberOfSim, int *resultArray)
{
	for (int simCount = 0; simCount<numberOfSim; simCount++)
	{
#if WIN32 | _WIN32
		float uniformRandom = float(rand()) / RAND_MAX;
#else
		float uniformRandom = drand48();
#endif
		int temp = (int)yFunction(getRandomVariable(uniformRandom));
		resultArray[simCount] =temp;
	}
}

double ExponentDistSimulation::getRandomVariable(double cdf)
{
	// CDF 0 means X is infinity
	return -log(cdf == 0 ? 0.0001f : cdf) / rateOfDistribution;//Not dividing by rate inside log as is provides same distribution scale as dividing outside
}

void ExponentDistSimulation::gatherNecessaryData()
{
	FunctionalSimulation::gatherNecessaryData();
	GET_STREAM_DATA("Enter Rate of Exponent : ", rateOfDistribution);
}

void BeerLambertLightPropagationSim::gatherNecessaryData()
{
	FunctionalSimulation::gatherNecessaryData();
	GET_STREAM_DATA("Rate of Light absorption per unit length of material simulating : ",rateOfAbsorbtion );
	GET_STREAM_DATA("Rate of Light scattering per unit length of material simulating : ", rateOfScatter);
	rateOfDistribution = rateOfAbsorbtion + rateOfScatter;
}



double PhysicalLightPropagationSim::getRandomVariable(double g)
{
#if WIN32 | _WIN32
	float uniformRandom = float(rand()) / RAND_MAX;
#else
	float uniformRandom = drand48();
#endif

	if (g != 0.0f)
	{
		double innerMost = (1 - g * g) / (1 - g + 2 * g*uniformRandom);
		return (1 + g * g - innerMost * innerMost) / (2 * g);
	}
	return 1 - 2 * uniformRandom;
}

PhysicalLightPropagationSim::~PhysicalLightPropagationSim()
{
	delete(randomLightPropagationDistance);
}

void PhysicalLightPropagationSim::gatherNecessaryData()
{
	nBounces = nAbsorbed = nTransmitted=0;
	if (randomLightPropagationDistance)
	{
		return;
	}

	randomLightPropagationDistance = new BeerLambertLightPropagationSim();
	randomLightPropagationDistance->gatherNecessaryData();
	
	pixelDensity = 1;
	slabSize = 10;
	slabHeight = 40;
	deflectionMag = 0.85f;
	photonGrpEnergy = 35;
	numberOfPasses = 1;
	rouletteSampleSize = 6;


	GET_STREAM_DATA("Pixel density of final image : ", pixelDensity);
	GET_TWO_STREAM_DATA("Enter Size and Thickness of slab : ", slabSize, slabHeight);
	GET_STREAM_DATA("Enter deflection intensity in range -1 to 1 : ", deflectionMag);
	GET_STREAM_DATA("Photon group initial energy : ", photonGrpEnergy);
	GET_STREAM_DATA("Number of Passes : ", numberOfPasses);
	GET_STREAM_DATA("Number of chances in Russian roulette : ", rouletteSampleSize);
	 

	imageSize = int(slabSize * pixelDensity);;
}

void PhysicalLightPropagationSim::fillSimulatedResult(int numberOfSim, int *resultArray)
{
	reflectedDataAllPass = new float[imageSize*imageSize];
	memset(reflectedDataAllPass, 0x0, sizeof(float) *imageSize*imageSize);

	transmittedDataAllPass = new float[imageSize*imageSize];
	memset(transmittedDataAllPass, 0x0, sizeof(float) *imageSize*imageSize);

	lightPropagatedData = new unsigned int[imageSize*imageSize];
	memset(lightPropagatedData, 0x0, sizeof(unsigned int) *imageSize*imageSize);
	lightReflectedData = new unsigned int[imageSize*imageSize];
	memset(lightReflectedData, 0x0, sizeof(unsigned int) *imageSize*imageSize);

	for (int i = 1; i <= numberOfPasses; i++)
	{
		nCurrentPass = i;
		doSimulationPass(numberOfSim);
		consoleOutResult(nullptr, 0);
	}

	delete[] reflectedDataAllPass;
	delete[] transmittedDataAllPass;
	delete[] lightPropagatedData;
	delete[] lightReflectedData;

	reflectedDataAllPass = transmittedDataAllPass = nullptr;
	lightPropagatedData = lightReflectedData = nullptr;
}

void PhysicalLightPropagationSim::doSimulationPass(int nPhotons)
{
	float slabExtend = slabSize / 2;

	for (int i = 0; i < nPhotons; i++)
	{
		Vector3 position = enteringPoint;
		Vector3 direction = enteringDirection;
		float energy = photonGrpEnergy;

		while (1)
		{
#if WIN32 | _WIN32
			float uniformRandom = float(rand()) / RAND_MAX;
#else
			float uniformRandom = drand48();
#endif
			double distance = randomLightPropagationDistance->getRandomVariable(uniformRandom);

			float s = (direction.getZ() > 0 ? slabHeight - position.getZ() : -position.getZ()) / direction.getZ();

			if (distance > s)// Photon is outside the slab
			{
				position = position + direction * s;
				float *fillPtr;
				unsigned int *imageDataPtr;
				if (direction.getZ() > 0)
				{
					fillPtr = transmittedDataAllPass;
					imageDataPtr = lightPropagatedData;
					nTransmitted++;
				}
				else
				{
					fillPtr = reflectedDataAllPass;
					imageDataPtr = lightReflectedData;
					nBounces++;
				}
				if (position.getX() >= -slabExtend && position.getY() >= -slabExtend && position.getX() <= slabExtend &&
					position.getY() <= slabExtend)
				{
					int pixelWidth = int((position.getX() + slabExtend)*pixelDensity);
					int pixelHeight = int((position.getY() + slabExtend)*pixelDensity);
					fillPtr[int(pixelHeight*imageSize) + pixelWidth] += energy / photonGrpEnergy;
					imageDataPtr[int(pixelHeight*imageSize) + pixelWidth] =
						(int)min(255.0f, (fillPtr[int(pixelHeight*imageSize) + pixelWidth] / nCurrentPass) * 255);
				}
				break;
			}
			position = position + direction * distance;

			float energyLoss = randomLightPropagationDistance->getAbsorptionRate() / randomLightPropagationDistance->getTransmitFallOffRate();
			energy -= energyLoss;
			energy = max(0.0f, energy);
			if (energy <= 0.01f)
			{
#if WIN32 | _WIN32
				uniformRandom = float(rand()) / RAND_MAX;
#else
				uniformRandom = drand48();
#endif
				if (uniformRandom >= 1.0f / rouletteSampleSize)
				{
					energy *= rouletteSampleSize;
				}
				else
				{
					nAbsorbed++;
					break;
				}
			}

			deflectPhoton(direction, min(max(deflectionMag, -1.0f), 1.0f));
		}
	}
}

void PhysicalLightPropagationSim::deflectPhoton(Vector3& direction, float g)
{
	double cosTheta = getRandomVariable(g);
	double sinTheta = sqrt(1 - cosTheta * cosTheta);

#if WIN32 | _WIN32
	float uniformRandom = float(rand()) / RAND_MAX;
#else
	float uniformRandom = drand48();
#endif

	float phi = float(uniformRandom * 2 * M_PI);
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);
	double x, y, z;

	if (abs(direction.getZ()) == 1.0f)
	{
		x = sinTheta * cosPhi;
		y = direction.getZ()*sinTheta* sinPhi;
		z = direction.getZ()* cosTheta;
	}
	else {
		double oneMinusZSquare = 1 - direction.getZ()*direction.getZ();
		double innerSoln = sinTheta / sqrt(oneMinusZSquare);

		x = cosTheta * direction.getX() + innerSoln * (cosPhi*direction.getX()*direction.getZ() - sinPhi * direction.getY());
		y = cosTheta * direction.getY() + innerSoln * (cosPhi*direction.getY()*direction.getZ() + sinPhi * direction.getX());
		z = cosTheta * direction.getZ() - innerSoln * cosPhi*oneMinusZSquare;
	}


	direction.setX((float)x);
	direction.setY((float)y);
	direction.setZ((float)z);
}

float* PhysicalLightPropagationSim::smoothImageData(const unsigned int *imgData, int width, int height, int smoothSize)
{
	assert(width > smoothSize && height > smoothSize);
	float* intermediateData = new float[width*height];
	memset(intermediateData, 0x0, sizeof(float)*width*height);

	float* returnData = new float[width*height];
	memset(returnData, 0x0, sizeof(float)*width*height);


	smoothSize = (smoothSize >> 1) << 1 == smoothSize ? smoothSize + 1 : smoothSize;
	int halfSmooth = smoothSize / 2;

	for (int h = 0; h < height; h++)
	{
		int hIdx = h * width;
		for (int w = 0; w < width; w++)
		{
			int wIdx = hIdx + w;
			int c = w - halfSmooth;
			if (c < 0)
			{
				int d = halfSmooth - w;
				intermediateData[wIdx] += d * imgData[hIdx];
				c = hIdx;
				while (c <= wIdx + halfSmooth)
				{
					intermediateData[wIdx] += imgData[c];
					c++;
				}
			}
			else if (w + halfSmooth >= width)
			{
				int d = (w + halfSmooth) + 1 - width;
				intermediateData[wIdx] += d * imgData[hIdx + width - 1];
				c += hIdx;
				while (c - hIdx < width && c <= wIdx + halfSmooth)
				{
					intermediateData[wIdx] += imgData[c];
					c++;
				}
			}
			else
			{
				c += hIdx;
				int subIdx = max(c - 1, hIdx);
				int addIdx = min(hIdx + width - 1, wIdx + halfSmooth);
				intermediateData[wIdx] = intermediateData[wIdx - 1] - imgData[subIdx] + imgData[addIdx];
			}
		}
	}

	for (int w = 0; w < width; w++)
	{
		for (int h = 0; h < height; h++)
		{
			int hIdx = h * width;
			int wIdx = hIdx + w;

			int c = h - halfSmooth;

			if (c < 0)
			{
				int d = halfSmooth - h;
				returnData[wIdx] += d * intermediateData[w];
				c = w;
				while (c <= wIdx + halfSmooth * width)
				{
					returnData[wIdx] += intermediateData[c];
					c += width;
				}
			}
			else if (h + halfSmooth >= height)
			{
				int d = (h + halfSmooth) - height + 1;
				returnData[wIdx] += d * intermediateData[w + (height - 1)*width];
				c = w + c * width;
				while ((c - w) / width  < height && c <= wIdx + halfSmooth * width)
				{
					returnData[wIdx] += intermediateData[c];
					c += width;
				}
			}
			else
			{
				c = w + c * width;
				int subIdx = max(c - width, w);
				int addIdx = min(w + (height - 1)*width, wIdx + halfSmooth * width);
				returnData[wIdx] = returnData[wIdx - width] - intermediateData[subIdx] + intermediateData[addIdx];
			}
		}
	}
	for (int h = 0; h < height; h++)
	{
		int hIdx = h * width;
		for (int w = 0; w < width; w++)
		{
			int wIdx = hIdx + w;
			returnData[wIdx] /= smoothSize * smoothSize;
		}
	}

	delete[] intermediateData;


	return returnData;
}

void PhysicalLightPropagationSim::consoleOutResult(int *result, int size)
{
	if (!lightPropagatedData || !lightReflectedData)
		return;
	printf("After Pass %d\n", nCurrentPass);

	string fileOuput = "C:/Users/JesJas/Documents/Visual Studio 2017/CppProjects/3dLearning/PlotData/lightPropagation%d.ppm";
	
	char *name=new char[fileOuput.length()+1];

	sprintf_s(name, fileOuput.length() + 1, fileOuput.c_str(), nCurrentPass);

	char *imageData = new char[imageSize*imageSize * 3];
	float* smoothedData = smoothImageData(lightPropagatedData, imageSize, imageSize, 10);
	for (int i = 0; i < imageSize*imageSize; i++)
	{
		imageData[i * 3] = imageData[i * 3 + 1] = imageData[i * 3 + 2] = (char)smoothedData[i];
	}
	delete[]smoothedData;

	if (!utils::writeToImage(name, imageSize, imageSize,imageData))
	{
		perror("Error in writing File");
		delete[]name;
		return;
	}
	delete[]name;

	// Reflectance writing
	fileOuput = "C:/Users/JesJas/Documents/Visual Studio 2017/CppProjects/3dLearning/PlotData/lightReflected%d.ppm";

	name = new char[fileOuput.length() + 1];

	sprintf_s(name, fileOuput.length() + 1, fileOuput.c_str(), nCurrentPass);

	smoothedData = smoothImageData(lightReflectedData, imageSize, imageSize, 10);
	for (int i = 0; i < imageSize*imageSize; i++)
	{
		imageData[i * 3] = imageData[i * 3 + 1] = imageData[i * 3 + 2] = (char)smoothedData[i];
	}
	delete[]smoothedData;

	if (!utils::writeToImage(name, imageSize, imageSize, imageData))
	{
		perror("Error in writing File");
		delete[]name;
		return;
	}

	delete[]name;

	printf("No of Bounced Photons : %d , Number absorbed : %d , Number Transmitted : %d\n", nBounces, nAbsorbed, nTransmitted);
}
