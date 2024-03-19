#include "MyVector.h"
/*
__forceinline MyVector::MyVector()
{
	vec[0] = 0;
	vec[1] = 0;
	vec[2] = 0;
}

__forceinline MyVector::MyVector(const MyVector& vec2)
{
	vec[0] = vec2.vec[0];
	vec[1] = vec2.vec[1];
	vec[2] = vec2.vec[2];
}

__forceinline MyVector::MyVector(const double vec2[3])
{
	vec[0] = vec[0];
	vec[1] = vec[1];
	vec[2] = vec[2];
}

__forceinline double MyVector::absValue()
{
	double sqrSum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	return pow(sqrSum, 0.5);
}

__forceinline MyVector MyVector::vecMul(MyVector vec2)
{
	MyVector res;
	res.vec[0] = vec[1] * vec2.vec[2] - vec[2] * vec2.vec[1];
	res.vec[1] = vec[2] * vec2.vec[0] - vec[0] * vec2.vec[2];
	res.vec[2] = vec[0] * vec2.vec[1] - vec[1] * vec2.vec[0];
	return res;
}

__forceinline MyVector& MyVector::operator=(const MyVector& vec2)
{
	vec[0] = vec2.vec[0];
	vec[1] = vec2.vec[1];
	vec[2] = vec2.vec[2];
	return (*this);
}

__forceinline double& MyVector::operator[](const int& num)
{
	return vec[num];
}

__forceinline const double& MyVector::operator[](const int& num) const
{
	return vec[num];
}

__forceinline MyVector MyVector::operator+(const MyVector& vec2)
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] += vec2[i];
	}
	return res;
}

__forceinline MyVector MyVector::operator*(const double& op2) const
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] *= op2;
	}
	return res;
}

__forceinline MyVector MyVector::operator/(const double& op2) const
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] /= op2;
	}
	return res;
}

//bool MyVector::operator==(const MyVector& vec2) const
//{
//	return *this == vec2;
//}
*/