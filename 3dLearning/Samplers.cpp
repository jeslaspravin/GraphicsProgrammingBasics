#include "stdafx.h"
#include "Samplers.h"
#include<cstdio>
#include<random>
#include<cstdlib>
#include <map>
#include <iterator>

void StatisticSampler::calculateMeanVariance()
{
	int elementryPopCnt,maxPopWeight,minPopWeight,numOfSamples;

	GET_STREAM_DATA("Enter population elementry count : ", elementryPopCnt);
	GET_TWO_STREAM_DATA("Enter min and max population per item : ", minPopWeight, maxPopWeight);
	GET_STREAM_DATA("Number of sampling needs to be done : ", numOfSamples);

	mt19937 rng;
	rng.seed(1117);

	int popSize=0;
	int *populationDistribution=new int[elementryPopCnt];

	double expectedValueMean=0, expectedVariance=0;

	uniform_int_distribution<uint32_t> popWeightGen(minPopWeight, maxPopWeight);

	for (int popItemIndex = 0; popItemIndex < elementryPopCnt; popItemIndex++)
	{
		int popWeight = popWeightGen(rng);
		populationDistribution[popItemIndex] = popWeight;
		popSize += popWeight;
		expectedValueMean += popWeight * popItemIndex;// Calculating Numerator part alone sum(value * weightage)
		expectedVariance += popWeight * popItemIndex*popItemIndex;// Calculating Numerator of first part of Variance eqn E(X^2)
	}

	expectedValueMean /= popSize;// Averaging all weights
	expectedVariance /= popSize;// Averaging for first part - E(X^2)
	expectedVariance -= expectedValueMean * expectedValueMean;// Second part E(X)^2


	uniform_int_distribution<uint32_t> sampleSizeGen(elementryPopCnt-10, elementryPopCnt - 1);//-1 since index goes -1 of size
	uniform_int_distribution<uint32_t> weightScaledOutcomeSelector(0, popSize - 1);

	double avgMean = 0, avgVar = 0;

	for (int samplerCount = 0; samplerCount < numOfSamples; samplerCount++)
	{
		double sampleMean = 0, sampleVariance = 0;
		int numToSample = sampleSizeGen(rng);
		for (int sampleCount = 0; sampleCount < numToSample; sampleCount++)
		{
			int selectedWeight = weightScaledOutcomeSelector(rng);
			int sampleSpaceIndex;
			for (sampleSpaceIndex = 0; sampleSpaceIndex < elementryPopCnt; sampleSpaceIndex++)
			{
				selectedWeight -= populationDistribution[sampleSpaceIndex];
				if (selectedWeight < 0)
					break;
			}
			sampleMean += sampleSpaceIndex;
			sampleVariance += sampleSpaceIndex * sampleSpaceIndex;
		}
		sampleMean /= numToSample;
		sampleVariance /= numToSample;
		sampleVariance -= sampleMean * sampleMean;

		printf("%d Sampler Pass :\n\tSample Mean : %0.5f, Sample Variance : %0.5f \n\t of %d Samples\n", samplerCount+1, sampleMean,
			sampleVariance, numToSample);

		avgMean = (avgMean + sampleMean) / 2;
		avgVar = (avgVar + sampleVariance) / 2;
	}

	printf("Population %d had Expected Mean %0.4f ,Variance %0.4f\n", popSize, expectedValueMean, expectedVariance);
	printf("Average Mean and Variance after %d Sampling : %0.5f,%0.5f", numOfSamples, avgMean, avgVar);

	delete[] populationDistribution;
}


void StatisticSampler::calculateMeanVarianceCount()
{
	int elementryPopCnt, maxPopWeight, minPopWeight, numOfSamples;

	GET_STREAM_DATA("Enter population elementry count : ", elementryPopCnt);
	GET_TWO_STREAM_DATA("Enter min and max population per item : ", minPopWeight, maxPopWeight);
	GET_STREAM_DATA("Number of sampling needs to be done : ", numOfSamples);

	mt19937 rng;
	rng.seed(1117);

	int popSize = 0;
	int *populationDistribution = new int[elementryPopCnt];

	double expectedValueMean = 0, expectedVariance = 0;

	uniform_int_distribution<uint32_t> popWeightGen(minPopWeight, maxPopWeight);

	for (int popItemIndex = 0; popItemIndex < elementryPopCnt; popItemIndex++)
	{
		int popWeight = popWeightGen(rng);
		populationDistribution[popItemIndex] = popWeight;
		popSize += popWeight;
		expectedValueMean += popWeight * popItemIndex;// Calculating Numerator part alone sum(value * weightage)
		expectedVariance += popWeight * popItemIndex*popItemIndex;// Calculating Numerator of first part of Variance eqn E(X^2)
	}

	expectedValueMean /= popSize;// Averaging all weights
	expectedVariance /= popSize;// Averaging for first part - E(X^2)
	expectedVariance -= expectedValueMean * expectedValueMean;// Second part E(X)^2


	uniform_int_distribution<uint32_t> weightScaledOutcomeSelector(0, popSize - 1);

	double avgMean = 0, avgVar = 0;
	int numToSample;
	GET_STREAM_DATA("Enter number of Samples : ", numToSample);

	map<int, int> eachCardData;
	for (int cardNo = 1; cardNo <= elementryPopCnt; ++cardNo)
	{
		eachCardData.insert(pair<int, int>(cardNo, 0));
	}

	for (int samplerCount = 0; samplerCount < numOfSamples; samplerCount++)
	{
		double sampleMean = 0, sampleVariance = 0;
		for (int sampleCount = 0; sampleCount < numToSample; sampleCount++)
		{
			int selectedWeight = weightScaledOutcomeSelector(rng);
			int sampleSpaceIndex;
			for (sampleSpaceIndex = 0; sampleSpaceIndex < elementryPopCnt; sampleSpaceIndex++)
			{
				selectedWeight -= populationDistribution[sampleSpaceIndex];
				if (selectedWeight < 0)
					break;
			}
			sampleMean += sampleSpaceIndex;
			sampleVariance += sampleSpaceIndex * sampleSpaceIndex;
		}
		sampleMean /= numToSample;
		sampleVariance /= numToSample;
		sampleVariance -= sampleMean * sampleMean;

		int roundedCardVal = (int)sampleMean + 1;
		eachCardData[roundedCardVal]++;

		printf("%d Sampler Pass :\n\tSample Mean : %0.5f, Sample Variance : %0.5f \n\t of %d Samples\n", samplerCount + 1, sampleMean,
			sampleVariance, numToSample);

		avgMean += sampleMean;
		avgVar += sampleMean* sampleMean;
	}

	avgMean /= numOfSamples;
	avgVar /= numOfSamples;
	avgVar -= avgMean * avgMean;

	printf("Population %d had Expected Mean %0.4f ,SD %0.4f\n", popSize, expectedValueMean, sqrt(expectedVariance));
	printf("Average Mean and SD after %d Sampling : %0.5f\n", numOfSamples, avgMean);
	for (int i = 0; i < 50; i++) cout << "-";
	cout << '\n';
	for (auto pairItem = eachCardData.begin(); pairItem != eachCardData.end(); ++pairItem)
	{
		cout << pairItem->first<<((pairItem->second>0)?"\t|[":"\t|");
		for (int i = 0; i < pairItem->second; i+=5)
		{
			cout << "|";
		}
		cout << ((pairItem->second > 0) ? "]":"")<<'\n';
	}
	for (int i = 0; i < 50; i++) cout << "-";
	delete[] populationDistribution;
}
