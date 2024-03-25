#include "MyVector.h"

MyVector::MyVector()
{
	vec[0] = 0;
	vec[1] = 0;
	vec[2] = 0;
}

MyVector::MyVector(const MyVector& vec2)
{
	vec[0] = vec2.vec[0];
	vec[1] = vec2.vec[1];
	vec[2] = vec2.vec[2];
}

MyVector::MyVector(const double vec2[3])
{
	vec[0] = vec[0];
	vec[1] = vec[1];
	vec[2] = vec[2];
}

double MyVector::absValue()
{
	double sqrSum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	return pow(sqrSum, 0.5);
}

MyVector MyVector::vecMul(MyVector vec2)
{
	MyVector res;
	res.vec[0] = vec[1] * vec2.vec[2] - vec[2] * vec2.vec[1];
	res.vec[1] = vec[2] * vec2.vec[0] - vec[0] * vec2.vec[2];
	res.vec[2] = vec[0] * vec2.vec[1] - vec[1] * vec2.vec[0];
	return res;
}

MyVector& MyVector::operator=(const MyVector& vec2)
{
	vec[0] = vec2.vec[0];
	vec[1] = vec2.vec[1];
	vec[2] = vec2.vec[2];
	return (*this);
}

double& MyVector::operator[](const int& num)
{
	return vec[num];
}

const double& MyVector::operator[](const int& num) const
{
	return vec[num];
}

MyVector MyVector::operator+(const MyVector& vec2)
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] += vec2[i];
	}
	return res;
}

MyVector MyVector::operator*(const double& op2) const
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] *= op2;
	}
	return res;
}

MyVector MyVector::operator/(const double& op2) const
{
	MyVector res(*this);
	for (int i = 0; i < 3; i++)
	{
		res[i] /= op2;
	}
	return res;
}

std::ostream& operator<<(std::ostream& os, const MyVector& vec)
{
	os << vec[0] << ' ' << vec[1] << ' ' << vec[2] << '\n';
	return os;
}
//bool MyVector::operator==(const MyVector& vec2) const
//{
//	return *this == vec2;
//}
