#include "stdafx.h"
#include "Integration.h"
#include <algorithm>

void colors::convertSPDtoXYZbyMC(const double *spd, const int spdSize, Vector3 &xyz, int nSamples)
{
	float brightnessNormalize = 0;
	xyz = Vector3();
	for (int i = 0; i < nSamples; i++)
	{
#if WIN32 | _WIN32
		float random = float(rand()) / RAND_MAX;
#else// #if WIN32 | _WIN32
		float random = drand48();
#endif// #if WIN32 | _WIN32

#if USE_STRATIFIED_SAMPLES
		random = (i + random) / nSamples;
#endif// #if USE_STRATIFIED_SAMPLES

		float lambda = (VISIBLE_SPECTRUM_MAX - VISIBLE_SPECTRUM_MIN)*random;
		float spdDataIdx = (spdSize - 1)*random;

		float spdValue = (float)utils::lerp(spd[(int)spdDataIdx], spd[min(int(spdDataIdx) + 1, spdSize - 1)], 
			spdDataIdx - (int)spdDataIdx);

		float colorMatchIdx = lambda / CIE_FN_LAMBDA_DELTA;
		int lowerIdx = (int)colorMatchIdx;
		int higherIdx = min(lowerIdx + 1, N_COLOR_BINS * 2 - 1);
		float dt = colorMatchIdx - lowerIdx;

		xyz += Vector3(
			float(spdValue*utils::lerp(COLOR_MATCHING_FN[0][lowerIdx], COLOR_MATCHING_FN[0][higherIdx], dt)),
			float(spdValue*utils::lerp(COLOR_MATCHING_FN[1][lowerIdx], COLOR_MATCHING_FN[1][higherIdx], dt)),
			float(spdValue*utils::lerp(COLOR_MATCHING_FN[2][lowerIdx], COLOR_MATCHING_FN[2][higherIdx], dt))
			);

		brightnessNormalize += (float)utils::lerp(COLOR_MATCHING_FN[1][lowerIdx], COLOR_MATCHING_FN[1][higherIdx], dt);
	}

	brightnessNormalize *= (VISIBLE_SPECTRUM_MAX - VISIBLE_SPECTRUM_MIN)/nSamples;
	xyz *= (VISIBLE_SPECTRUM_MAX - VISIBLE_SPECTRUM_MIN) / (nSamples * brightnessNormalize);
}

void colors::createMcBertChart(int tileSize, int nSamples)
{
	static const float gamma = 1/2.2f ;
	int imgWidth = tileSize * 6,imgHeight= tileSize * 4;
	char* imageData = new char[imgWidth*imgHeight * 3];
	float* imageBufferData = new float[imgWidth*imgHeight * 3];
	memset(imageBufferData, 0x0, sizeof(float)*imgWidth*imgHeight * 3);

	string imagePath = "C:/Users/JesJas/Documents/Visual Studio 2017/CppProjects/3dLearning/PlotData/McBerthChart.ppm";
	string passImagePath = "C:/Users/JesJas/Documents/Visual Studio 2017/CppProjects/3dLearning/PlotData/McBerthChartPass%d.ppm";
	const char* imgPathPtr = imagePath.c_str();

	int nPass;
	GET_STREAM_DATA("Enter Pass Count : ", nPass);
	for(int i = 1 ;i<=nPass;i++)
	{
		for (int h = 0; h < imgHeight; h++)
		{
			int macBerthRow = int(h / tileSize);
			int processedIndex = h * imgWidth;
			for (int w = 0; w < imgWidth; w++)
			{
				int macBerthColumn = int(w / tileSize);
				int macBerthSpdIdx = macBerthRow * 6 + macBerthColumn;
				int imgIdx = processedIndex + w;
				imgIdx *= 3;

				Vector3 xyz, rgb;
				convertSPDtoXYZbyMC(MCBERTH_SPECTRAL_DATA[macBerthSpdIdx], N_COLOR_BINS, xyz, nSamples);
				convertXYZtoRGB(xyz, rgb);
				imageBufferData[imgIdx] += rgb.get(0);
				imageBufferData[imgIdx + 1] += rgb.get(1);
				imageBufferData[imgIdx + 2] += rgb.get(2);
				// Gamma corrected sRGB is stored in file format
				imageData[imgIdx] =char( powf(min(1.f, imageBufferData[imgIdx] / i), gamma)*255.f);
				imageData[imgIdx + 1] = char(powf(min(1.f, imageBufferData[imgIdx+1] / i),gamma)*255.f);
				imageData[imgIdx + 2] = char(powf(min(1.f, imageBufferData[imgIdx+2] / i), gamma)*255.f);
			}
		}
		if (i < nPass)
		{
			char* passPath=new char[passImagePath.length() + 1];;
			sprintf_s(passPath, passImagePath.length() + 1, passImagePath.c_str(), i);
			utils::writeToImage(passPath, imgWidth, imgHeight, imageData);
		}
	}

	utils::writeToImage(imgPathPtr, imgWidth, imgHeight, imageData);
	delete[] imageData;
}
