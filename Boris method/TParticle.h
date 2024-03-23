#pragma once
#include <iostream>
#include <vector>
#include "MyVector.h"
#include <fstream>

// CGS units
#define c 2.99792458e10

class TParticle
{
public:
	MyVector r_old; // Particle position (x, y, z)
	MyVector p_old; // Particle momentum (px, py, pz)
public:
	TParticle();
	TParticle(double _x, double _y, double _z,
		double _px, double _py, double _pz);

	//MyVector makeOneStep(const MyVector& E, const MyVector& B);
	MyVector getP();

};

