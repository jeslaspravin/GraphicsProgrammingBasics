#pragma once


struct Vector3
{
private:
	float x;
	float y;
	float z;

public:

	Vector3() :x(0), y(0), z(0) {}

	Vector3(float X,float Y,float Z):x(X),y(Y),z(Z)
	{}

	Vector3(float X,float Y):x(X),y(Y),z(0){}

	Vector3 operator *(Vector3 v);
	void operator *=(Vector3 v);
	Vector3 operator *(float value);
	void operator *=(float value);
	Vector3 operator *(double value);
	void operator *=(double value);
	Vector3 operator *(int value);
	void operator *=(int value);

	Vector3 operator +(Vector3 v);
	void operator +=(Vector3 v);

	void normalize();
	Vector3 getNormalized();
	float magnitude();

	float getX() const { return x; }
	void setX(float value) { x = value; }
	float getY() const { return y; }
	void setY(float value) { y = value; }
	float getZ() const { return z; }
	void setZ(float value) { z = value; }

	float get(int i)
	{
		switch (i)
		{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			perror("Invalid Index");
		}

		return -1;
	}

	void set(int i,float value)
	{
		switch (i)
		{
		case 0:
			x = value;
			break;
		case 1:
			y = value;
			break;
		case 2:
			z = value;
			break;
		default:
			perror("Invalid Index");
		}
	}
};
