// 3dLearning.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Distributions.h"
#include "Samplers.h"
#include "Simulations.h"
#include "Integration.h"
#include <cstdio>
#include <random>
#include <cstdlib>
#include <cmath>
#include <map>

void distribution()
{
	int opt = 0;
	std::cout << "Choose functionality:\n1.Binomial Distribution\n2.CDF for Normal Distribution\n0.Go to Main Menu\nYour selection : ";
	std::cin >> opt;
	Distributions* distribution = nullptr;
	switch (opt)
	{
	case 1:
		distribution = new BinomialDistribution();
		break;
	case 2:
		distribution = new NormalDistributionCDF();
		break;
	default:
		std::cout << "No valid option is selected\n";
	}

	if (distribution)
	{
		std::cout << "Probability of the observation is : " << distribution->calculateDistribution() << '\n';
		//distribution->printDistributionVsTrialCnt();
		delete distribution;
	}
}

void sampling()
{
	int opt = 0;
	std::cout << "Choose functionality:\n1.Random Sampler for Population\n0.Go to Main Menu\nYour selection : ";
	std::cin >> opt;

	BaseSampler* sampler = nullptr;

	switch (opt)
	{
	case 1:
		sampler = new StatisticSampler();
		break;
	default:
		std::cout << "No valid option is selected\n";
	}
	if (sampler)
	{
		sampler->calculateMeanVarianceCount();
		delete sampler;
	}
}

void simulating()
{
	int opt = 0;

	std::cout << "Choose functionality:\n1.Simulation with events around middle\n2.Exponential Simulation\n3.Beer-Lambert Light Propagation model\n4.Physical Light Propagation\n0.Go to Main Menu\nYour selection : ";
	std::cin >> opt;

	Simulations *sim = nullptr;

	switch (opt)
	{
	case 1:
		sim = new NormalDistributedArbitrarySim();
		break;
	case 2:
		sim = new ExponentDistSimulation();
		break;
	case 3:
		sim = new BeerLambertLightPropagationSim();
		break;
	case 4:
		sim = new PhysicalLightPropagationSim();
		break;
	default:
		std::cout << "No valid option is selected\n";
	}

	if (sim)
	{
		sim->gatherNecessaryData();
		sim->generateRequiredDistribution();
		int simCount = 0;
		GET_STREAM_DATA("Enter number of simulation to be done : ", simCount);
		int* simulationResults = new int[simCount];
		sim->fillSimulatedResult(simCount, simulationResults);
		sim->consoleOutResult(simulationResults,simCount);
		delete sim;
		delete[] simulationResults;
	}

}

void circleAreaHitMiss()
{
	float circleRadius;
	GET_STREAM_DATA("Enter Circle Radius : ", circleRadius);
	int numberOfSamples;
	GET_STREAM_DATA("Enter number of Samples to take : ", numberOfSamples);

	std::default_random_engine gen;
	std::uniform_real_distribution<float> distrb;
	
	int hits=0;

	for (int i = 0; i < numberOfSamples; i++)
	{
		float x = distrb(gen)*circleRadius;
		float y = distrb(gen)*circleRadius;
		float l = sqrt(x*x + y * y);
		if (l <= circleRadius)
			hits++;
	}

	printf("Area from Hit and Miss %0.3f , Area from formula %0.3f", circleRadius*circleRadius * 4*hits / numberOfSamples, M_PI*circleRadius*circleRadius);
}

void mcBerthChart()
{
	int tileSize, nSamples;
	GET_STREAM_DATA("Enter Tile Size in pixels : ", tileSize);
	GET_STREAM_DATA("Enter number of samples : ", nSamples);

	colors::createMcBertChart(tileSize, nSamples);
}

// Matrix inverse by Gauss-Jordan row elimination method
void invertMatrix();

template<typename BitType>
void printBits(BitType bits);

template<typename BitType>
void printBits(BitType bits)
{
	std::cerr << "Unsupported implementation";
}

template <>
void printBits<uint8_t>(uint8_t bitNum)
{
	for (uint8_t i = 1, mask = 1 << 7; i <= 8; i++)
	{
		std::cout << ((bitNum & mask) > 0) ? "1" : "0";
		mask = mask >> 1;
	}
	std::cout << std::endl;
}

template <>
void printBits<uint16_t>(uint16_t bitNum)
{
	for (uint16_t i = 1, mask = 1 << 15; i <= 16; i++)
	{
		std::cout << ((bitNum & mask) > 0) ? "1" : "0";
		mask = mask >> 1;
	}
	std::cout << std::endl;
}

template <>
void printBits<uint32_t>(uint32_t bitNum)
{
	for (uint32_t i = 1, mask = 1 << 31; i <= 32; i++)
	{
		std::cout << ((bitNum & mask) > 0) ? "1" : "0";
		mask = mask >> 1;
	}
	std::cout << std::endl;
}

uint32_t ieee754(float f);
float reverseIeee754(uint32_t binaryRep);

int main()
{

#if WIN32 | _WIN32
	srand(1117);
#else
	srand48(1117);
#endif

	int opt = 0;
	do {
		std::cout << "Choose Task:\n1.Distribution\n2.Sampling\n3.Simulation\n4.Circle Area by Hit and Miss Method\n5.McBerth Chart\n6.Matrix inverse by Gauss-Jordan row elimination method\n7.IEEE 754 float representation\n0.Exit\nYour selection : ";
		std::cin >> opt;

		switch (opt)
		{
		case 1:
			distribution();
			break;
		case 2:
			sampling();
			break;
		case 3:
			simulating();
			break;
		case 4:
			circleAreaHitMiss();
			break;
		case 5:
			mcBerthChart();
			break;
		case 6:
			invertMatrix();
			break;
		case 7:
		{
			float f;
			GET_STREAM_DATA("Enter float value to see IEEE 754 representation : ", f);
			std::cout << std::endl;
			uint32_t binaryRep = ieee754(f);
			std::cout << "\nIEEE 754 Representation ";
			printBits(binaryRep);
			std::cout << std::endl;
			std::cout << "\nDecoded value " << reverseIeee754(binaryRep) << std::endl << std::endl;
			break;
		}
		default:
			std::cout << "No valid option is selected,Quiting!\n";
		}
	} while (opt);
	
	/*
	 * Sample mean of subset sample space after large enough samples averaged starts settling to a value - Law of Large Number
	 */
	//std::mt19937 rng;
	//rng.seed(1117);

	//std::uniform_int_distribution<uint32_t> coin(0, 1);
	//int sum = 0;

	//for (int i = 0; i < 1000; i += 1)
	//{
	//	int toss = (int)coin(rng);
	//	sum += toss;
	//	printf("Tossed %d , %d : %0.2f\n", toss,i, float(sum) / i);
	//}
	return 0;
}


void invertMatrix()
{
	// Creating Matrix
	Matrix3x3f mat, m;
	for (uint8_t i = 0; i < 3; i++)
	{
		for (uint8_t j = 0; j < 3; j++)
		{
			float val;
			char buffer[] = "Enter element [%d][%d] of matrix : ";
			sprintf_s(buffer, buffer, i+1, j+1);
			GET_STREAM_DATA(buffer, val);
			m[i][j] = val;
		}
	}

	cout << "Input Matrix is :\n";
	for (uint8_t row = 0; row < 3; row++)
	{
		cout << "|\t";
		for (uint8_t column = 0; column < 3; column++)
		{
			cout << m[row][column] << '\t';
		}
		cout << "|\n";
	}

	// Start converting m to Identity matrix from column 0
	for (uint8_t column = 0; column < 3; column++)
	{
		// If Pivot of current column is zero swap pivot row with row with max absolute value in that column
		if (m[column][column] == 0)
		{
			uint8_t bigRowIdx=column;
			for (uint8_t row = 0; row < 3; row++)
			{
				bigRowIdx = m[row][column] != 0 && abs(m[row][column] > m[bigRowIdx][column]) ? column : bigRowIdx;
			}

			if (bigRowIdx == column)
			{
				perror("This Matrix cannot be inverted");
				return;
			}
			else
			{
				for (uint8_t c = 0; c < 3; c++)
				{
					swap(m[column][c], m[bigRowIdx][c]);
					swap(mat[column][c], mat[bigRowIdx][c]);
				}
			}
		}

		// Set each row of column to 0 except pivot
		for (uint8_t row = 0; row < 3; row++)
		{
			// Each element in the current row has to be subtracted by element in current column multiplied by corresponding element in 
			// row with index equal to column's index,current column element has to be divided by pivot of current column before multiplying.
			if (row != column)
			{
				float coeff = m[row][column] / m[column][column];
				for (uint8_t c = 0; c < 3; c++)
				{
					m[row][c] -= coeff * m[column][c];
					mat[row][c] -= coeff * mat[column][c];
				}
				// Ensure it is zero to avoid floating point precision
				m[row][column] = 0;
			}
		}
	}

	for (uint8_t row = 0; row < 3; row++)
	{
		float pivot = m[row][row];
		for (uint8_t column = 0; column < 3; column++)
		{
			// Now convert all diagonal to 1,Since we only need mat we can skip doing that to actual matrix m
			// m[row][column] /= pivot;
			mat[row][column] /= pivot;
		}
	}

	cout << "Inverted Matrix is :\n";
	for (uint8_t row = 0; row < 3; row++)
	{
		cout << "|\t";
		for (uint8_t column = 0; column < 3; column++)
		{
			cout << mat[row][column] << '\t\t';
		}
		cout << "|\n";
	}
}

uint32_t ieee754(float f)
{
	uint32_t isNegative = f < 0 ? 1 : 0;
	uint32_t binaryFormat = isNegative << 31;
	uint32_t intPart = (uint32_t)abs(f);
	float floatPart = abs(f) - intPart;

	uint32_t decimalBits = 0;
	uint32_t intBits = 0;
	int shiftCount = 0;


	// Finding 23 bit of float representation(mantissa)
	float fltPartTemp = floatPart;
	int bitPos = 31;
	while (fltPartTemp > std::numeric_limits<float>::epsilon() || bitPos >= 0)
	{
		float tmp = fltPartTemp * 2;
		uint32_t bit = (uint32_t)tmp;
		decimalBits |= (bit > 0) ? (1 << bitPos) : 0;
		bitPos--;
		fltPartTemp = tmp - bit;
	}

	std::cout << "Decimal Bits ";
	printBits(decimalBits);

	bitPos = 0;
	uint32_t tempPart = intPart;
	while (tempPart > 0)
	{
		uint32_t q = tempPart / 2;
		uint32_t r = tempPart % 2;
		intBits |= r << bitPos;
		bitPos++;
		tempPart = q;
	}

	std::cout << "Int Bits ";
	printBits(intBits);

	tempPart = intBits;
	while (((tempPart | 1) ^ 1) > 0)
	{
		uint32_t bit = tempPart & 1;
		decimalBits = (decimalBits >> 1) | (bit << 31);
		shiftCount++;
		tempPart = tempPart >> 1;
	}

	// Shifting to place properly in final float format
	decimalBits = decimalBits >> 9;

	std::cout << "Decimal Bits after shift ";
	printBits(decimalBits);

	// Middle 8bit of exponent value,127 is added as it is considered middle value representing 0 shift
	uint32_t exponent = 127 + shiftCount;

	// Shifting to place properly in final float format
	binaryFormat |= (exponent << 23) | decimalBits;

	return binaryFormat;
}

// Switch to use bit shift for decoding the float value
#define USE_BIT_SHIFT_ONLY 1
float reverseIeee754(uint32_t binaryRep)
{
	uint32_t mask9thBit = (1 << 8);
	// Shifting to get exponent that will be from 24 to 31 bits
	uint32_t exponentBits = binaryRep >> 23;
	uint8_t signBit = (mask9thBit & exponentBits) >> 8;
	exponentBits = exponentBits & ~(mask9thBit);
	std::cout << "Sign bit ";
	printBits(signBit);
	std::cout << "Exponent bits ";
	printBits(exponentBits);

	uint8_t exponent = (uint8_t)(exponentBits - 127);
	std::cout << "Exponent " << (int)exponent << std::endl;

	uint32_t decimalBits = binaryRep << 9;
	std::cout << "Decimal bits ";
	printBits(decimalBits);

	uint32_t mask = 1 << 31;
	// Decoding float can be done either by using 2^exponent and multiplying it to decimal from decimalBits
	// Or can be done only using bit shifts
#if USE_BIT_SHIFT_ONLY
	uint32_t intPart = 1;
	uint32_t decimalPart = decimalBits;
	for (uint8_t i = 0; i < exponent; i++)
	{
		intPart = (intPart << 1) | ((mask & decimalPart) == 0 ? 0 : 1);
		decimalPart = (decimalPart << 1);
	}
	std::cout << "Decimal bits after re shifting ";
	printBits(decimalPart);
	std::cout << "Whole number bits after re shifting ";
	printBits(intPart);

	float decodedValue = 0;

	// Whole number part
	for (uint8_t i = 0; intPart > 0; i++,intPart=(intPart>>1))
	{
		//decodedValue += (intPart & 1)*(float)pow(2, i);
		decodedValue += (intPart & 1)*(1<<i);
	}

	// fractional part
	for (uint8_t i = 1; decimalPart > 0; i++, decimalPart = (decimalPart << 1))
	{
		//decodedValue+= ((mask & decimalPart) == 0 ? 0 : 1)/(float)pow(2, i);
		decodedValue += ((mask & decimalPart) == 0 ? 0 : 1) / (float)(1<<i);
	}

	return (signBit ? -1 : 1)*decodedValue;
#else
	uint32_t tempBits = decimalBits;
	// Adding 1 to result here itself
	float intermediate = 1;
	for (int i = 1; tempBits != 0; i++)
	{
		intermediate += 1.0f / (float)pow(2, i) * ((mask & tempBits) == 0 ? 0 : 1);
		tempBits = tempBits << 1;
	}
	std::cout << "Intermediate multiplier " << intermediate << std::endl;

	return (signBit?-1:1)*intermediate * (float)pow(2, exponent);
#endif
}
