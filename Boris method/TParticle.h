#pragma once
#include <iostream>
#include <vector>
#include "MyVector.h"
#include <fstream>

#define c 299792458

class TParticle
{
	//Необходимость использование double или float требует дальнейшего исследования
	double m = 9e-28;
	double gamma_old;
	double gamma_new;
	double q = -4.8e-10;
	double delta_t = 1e-12;
	MyVector r_old; // Particle position (x, y, z)
	MyVector r_new;
	MyVector v_new;
	MyVector t;
	MyVector p_old; // Particle momentum (px, py, pz)
	MyVector pMinus;
	MyVector pPlus;
	MyVector p_new;
	MyVector s;
	MyVector pDeriv;
public:
	TParticle();
	TParticle(double _x, double _y, double _z,
		double _px, double _py, double _pz);

	MyVector makeOneStep(const MyVector& E, const MyVector& B);

};

