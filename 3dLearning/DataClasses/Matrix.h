#pragma once

template<typename T>
class Matrix4x4
{
private:
	T m[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

public:

	Matrix4x4(){}

	const T* operator [](uint8_t i) const
	{
		return m[i];
	}

	T* operator [] (uint8_t i)
	{
		return m[i];
	}
};

template<typename T>
class Matrix3x3
{
private:
	T m[3][3] = { { 1,0,0 },{ 0,1,0 },{ 0,0,1 }};

public:

	Matrix3x3() {}

	const T* operator [](uint8_t i) const
	{
		return m[i];
	}

	T* operator [] (uint8_t i)
	{
		return m[i];
	}
};


typedef Matrix4x4<float> Matrix4x4f;
typedef Matrix3x3<float> Matrix3x3f;