#pragma once
#include <iostream>
#include <math.h>

class MyVector
{
	double vec[3];
public:
	MyVector();
	MyVector(const MyVector& vec2);
	MyVector(const double vec2[3]);
	double absValue();
	MyVector vecMul(MyVector vec2);
	MyVector& operator=(const MyVector& vec2);
	double& operator[](const int& num);
	const double& operator[](const int& num) const;
	MyVector operator+(const MyVector& vec2);
	MyVector operator*(const double& op2) const;
	MyVector operator/(const double& op2) const;
	friend std::ostream& operator<<(std::ostream& os, const MyVector& vec);
	//bool operator==(const MyVector& vec2) const;
};

