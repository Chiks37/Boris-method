#pragma once
#include <iostream>
#include <math.h>

class MyVector
{
	double vec[3];
public:
	MyVector() {
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = 0;
	}
	MyVector(const MyVector& vec2) {
		vec[0] = vec2.vec[0];
		vec[1] = vec2.vec[1];
		vec[2] = vec2.vec[2];
	}
	MyVector(const double vec2[3]) {
		vec[0] = vec2[0];
		vec[1] = vec2[1];
		vec[2] = vec2[2];
	}
	double absValue() {
		double sqrSum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
		return sqrt(sqrSum);
	}

	MyVector vecMul(MyVector vec2) {
		MyVector res;
		res.vec[0] = vec[1] * vec2.vec[2] - vec[2] * vec2.vec[1];
		res.vec[1] = vec[2] * vec2.vec[0] - vec[0] * vec2.vec[2];
		res.vec[2] = vec[0] * vec2.vec[1] - vec[1] * vec2.vec[0];
		return res;
	}
	MyVector& operator=(const MyVector& vec2) {
		vec[0] = vec2.vec[0];
		vec[1] = vec2.vec[1];
		vec[2] = vec2.vec[2];
		return (*this);
	}
	double& operator[](const int& num) {
		return vec[num];
	}
	const double& operator[](const int& num) const {
		return vec[num];
	}
	MyVector operator+(const MyVector& vec2) {
		MyVector res(*this);
		res[0] += vec2[0];
		res[1] += vec2[1];
		res[2] += vec2[2];
		return res;
	}
	MyVector operator*(const double& op2) const {
		MyVector res(*this);
		res[0] *= op2;
		res[1] *= op2;
		res[2] *= op2;
		return res;
	}

	MyVector operator/(const double& op2) const {
		MyVector res(*this);
		res[0] /= op2;
		res[1] /= op2;
		res[2] /= op2;
		return res;
	}
	//bool operator==(const MyVector& vec2) const;
	friend std::ostream& operator<<(std::ostream& os, const MyVector& vec);

};

