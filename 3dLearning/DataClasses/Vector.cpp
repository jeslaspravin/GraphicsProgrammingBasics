#include "stdafx.h"

Vector3 Vector3::operator*(Vector3 v)
{
	return Vector3(x * v.x, y * v.y, z* v.z);
}

Vector3 Vector3::operator*(int value)
{
	return Vector3(float(x * value), float(y * value), float(z* value));
}

Vector3 Vector3::operator*(double value)
{
	return Vector3(float(x * value), float(y * value), float(z* value));
}

Vector3 Vector3::operator*(float value)
{
	return Vector3(x * value, y * value, z* value);
}

void Vector3::operator*=(Vector3 v)
{
	x *= v.x; y *= v.y; z *= v.z;
}

Vector3 Vector3::operator+(Vector3 v)
{
	return Vector3(x + v.x, y + v.y, z + v.z);
}

void Vector3::operator+=(Vector3 v)
{
	x += v.x; y += v.y; z += v.z;
}

void Vector3::operator*=(int value)
{
	x *= value; y *= value; z *= value;
}

void Vector3::operator*=(float value)
{
	x *= value; y *= value; z *= value;
}

void Vector3::operator*=(double value)
{
	x *= (float)value; y *= (float)value; z *= (float)value;
}

void Vector3::normalize()
{
	float mag = magnitude();
	x /= mag;
	y /= mag;
	z /= mag;
}

Vector3 Vector3::getNormalized()
{
	float mag = magnitude();
	return Vector3(x / mag, y / mag, z / mag);
}

float Vector3::magnitude()
{
	return sqrt(x*x + y * y + z * z);
}
