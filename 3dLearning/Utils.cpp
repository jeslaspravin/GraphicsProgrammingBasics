#include "stdafx.h"


bool const utils::writeToImage(const char* path ,int width, int height, char* data)
{

	FILE *file;
	fopen_s(&file, path, "wb");

	if (file == nullptr)
	{
		perror("Cannot open file");
		return false;
	}

	fprintf(file, "P6\n");
	fprintf(file, "%d %d\n", width, height);
	fprintf(file, "255\n");
	int writtenSize = fwrite(data, 1, height*width * 3, file);
	fclose(file);
	file = nullptr;
	return writtenSize == height * width * 3;
}

double utils::lerp(double minValue, double maxValue, float dt)
{
	return minValue + (maxValue - minValue)*dt;
}

float utils::lerp(float minValue, float maxValue, float dt)
{
	return minValue + (maxValue - minValue)*dt;
}
