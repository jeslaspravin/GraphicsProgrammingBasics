#include "stdafx.h"
#include "Distributions.h"
#include<cmath>

double BinomialDistribution::calculateDistribution()
{
	float outComeProbability = 0;
	int numberOfSampling = 1, expectedCount = 1;
	cout << "Enter probability of an outcome of a Binomial Sampling : ";
	cin >> outComeProbability;
	cout << "Enter number of sampling done : ";
	cin >> numberOfSampling;
	cout << "Enter expected number of observation of current item : ";
	cin >> expectedCount;

	double cN = Utils::fact(numberOfSampling) / (float)(Utils::fact(expectedCount)*Utils::fact(numberOfSampling - expectedCount));
	return cN * pow(outComeProbability, expectedCount)*pow(1 - outComeProbability, numberOfSampling - expectedCount);
}

void BinomialDistribution::printDistributionVsTrialCnt()
{
	float outComeProbability = 0;
	int minNumberOfTrials, maxNumberOfTrials,stepSize;
	float ratioOfOutcomeCnt;
	cout << "Enter probability of an outcome of a Binomial Sampling : ";
	cin >> outComeProbability;
	cout << "Enter Min and Max Trails count : ";
	cin >> minNumberOfTrials>>maxNumberOfTrials;
	cout << "Enter Step size to increment Trails : ";
	cin >> stepSize; 
	cout << "Enter ratio of Outcome count from current trails count : ";
	cin >> ratioOfOutcomeCnt;

	for (int trailCnt = minNumberOfTrials; trailCnt <= maxNumberOfTrials; trailCnt += stepSize)
	{
		int expectedCount = (int)(ratioOfOutcomeCnt*trailCnt);
		//int expectedCount = (int)(ratioOfOutcomeCnt);
		double cN = Utils::fact(trailCnt) / (double)(Utils::fact(expectedCount)*Utils::fact(trailCnt - expectedCount));
		double probability=cN * pow(outComeProbability, expectedCount)*pow(1 - outComeProbability, trailCnt - expectedCount);
		cout << "Probability when trail count " << trailCnt << " : " << probability << '\n';
	}
}

inline double Utils::fact(int factNumber)
{
	return factNumber > 0 ? factNumber* fact(factNumber-1): 1;
}

double NormalDistributionCDF::calculateDistribution()
{
	float mean, sd,cdfParamValue;
	GET_TWO_STREAM_DATA("Enter Mean and Standard deviation of Normal Distribution : ", mean, sd);
	
	int minBound, maxBound, numOfSamples;
	GET_TWO_STREAM_DATA("Enter Minimum and Maximum bound of parameter range of population : ", minBound, maxBound);
	GET_STREAM_DATA("Enter number of samples to take from range : ", numOfSamples);
	GET_STREAM_DATA("Enter parameter value for which CDF is needed : ", cdfParamValue);

	double dx = (maxBound - minBound) / (double)numOfSamples;

	double cdf = 0;
	float currParamValue = (float)minBound;
	for (int i = 0; i < numOfSamples && currParamValue<cdfParamValue; i++)
	{
		// Using Mid Raimann Sum Method to find integral between min and max range
		currParamValue = minBound + (maxBound-minBound)*(i + 0.5f) / numOfSamples;
		float currClampedParamValue = currParamValue > cdfParamValue ? cdfParamValue : currParamValue;

		double pr = exp(-1 * (currClampedParamValue - mean)*(currClampedParamValue - mean) / (2 * sd*sd)) / (sd * sqrt(2 * M_PI));
		cdf += pr*dx;

		printf("CDF for Parameter Value %0.2f is %0.5f\n", currClampedParamValue, cdf);
	}
	return cdf;
}

void NormalDistributionCDF::printDistributionVsTrialCnt()
{
	float mean, sd;
	GET_TWO_STREAM_DATA("Enter Mean and Standard deviation of Normal Distribution : ", mean, sd);

	int minBound, maxBound, numOfSamples;
	GET_TWO_STREAM_DATA("Enter Minimum and Maximum bound of parameter range of population : ", minBound, maxBound);
	GET_STREAM_DATA("Enter number of samples to take from range : ", numOfSamples);

	double dx = (maxBound - minBound) / (double)numOfSamples;

	double cdf = 0;
	for (int i = 0; i < numOfSamples; i++)
	{
		// Using Mid Raimann Sum Method to find integral between min and max range
		float currParamValue = minBound + (maxBound - minBound)*(i + 0.5f) / numOfSamples;

		double pr = exp(-1*(currParamValue - mean)*(currParamValue - mean) / (2 * sd*sd)) / (sd * sqrt(2 * M_PI));
		cdf += pr*dx;

		printf("CDF for Parameter Value %0.2f is %0.5f\n", currParamValue, cdf);
	}
	printf("Final CDF is %0.5f\n", cdf);
}
