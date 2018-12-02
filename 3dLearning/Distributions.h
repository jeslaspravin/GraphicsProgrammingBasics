#pragma once


class Distributions
{
public:
	virtual double calculateDistribution() = 0;

	virtual void printDistributionVsTrialCnt() = 0;
};

class BinomialDistribution : public Distributions
{
public:
	virtual double calculateDistribution() override;
	virtual void printDistributionVsTrialCnt() override;
};

class NormalDistributionCDF : public Distributions
{
public:
	virtual double calculateDistribution() override;
	virtual void printDistributionVsTrialCnt() override;
};


namespace Utils
{
	inline double fact(int factNumber);
}