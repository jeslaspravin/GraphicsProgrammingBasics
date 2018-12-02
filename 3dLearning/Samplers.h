#pragma once

class BaseSampler
{
public:
	virtual void calculateMeanVariance() = 0;
	virtual void calculateMeanVarianceCount() = 0;
};

class StatisticSampler :public BaseSampler
{
	virtual void calculateMeanVariance() override;
	virtual void calculateMeanVarianceCount() override;
};