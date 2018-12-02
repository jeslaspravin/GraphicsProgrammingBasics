#pragma once


namespace utils {

	bool const writeToImage(const char* path, int width, int height, char* data);


	inline bool const writeToImage(const char* path, int width, int height, int* data)
	{
		char* imageCharData = new char[width*height * 3];
		for (int h = 0; h < height; h++)
		{
			for (int w = 0; w < width; w++)
			{
				int index = h * width + w;
				imageCharData[index] = (char)data[index];
			}
		}
		bool bResult = utils::writeToImage(path, width, height, imageCharData);

		delete []imageCharData;

		return bResult;
	}

	double lerp(double minValue, double maxValue, float dt);

	float lerp(float minValue, float maxValue, float dt);

}