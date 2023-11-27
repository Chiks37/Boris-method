#pragma once
#include <iostream>
#include <vector>
class MyVector
{
	std::vector<double> vec;
public:
	MyVector();
	MyVector(const MyVector& vec2);
	double absValue();
	MyVector vecMul(MyVector vec2);
	MyVector& operator=(std::vector<double> vec2);
	double& operator[](const int& num);
	const double& operator[](const int& num) const;
	MyVector operator+(const MyVector& vec2);
	MyVector operator*(const double& op2) const;
	MyVector operator/(const double& op2) const;
};
