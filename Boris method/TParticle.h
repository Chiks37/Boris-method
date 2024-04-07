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
	//Необходимость использование double или float требует дальнейшего исследования
	int partsCount;
	//double m = 9.10958215e-28;  // use more exact constants
	//double* gamma_old;
	//double* gamma_new;
	//double q = -4.803e-10;
	//double delta_t = 1e-12;
	MyVector* r_old; // Particle position (x, y, z)
	//MyVector* r_new;
	//MyVector* v_old;
	//MyVector* v_new;
	//MyVector* t;
	MyVector* p_old; // Particle momentum (px, py, pz)
	//MyVector* pMinus;
	//MyVector* pPlus;
	//MyVector* p_new;
	//MyVector* s;
	//MyVector* pDeriv;

	//double mcbased = m * c;
	//double qdthbased = q * delta_t * 0.5;

	TParticle(double _x, double _y, double _z,
		double _px, double _py, double _pz, int _partsCount);
	TParticle(int _partsCount);
	TParticle(const TParticle& ob2);
	//MyVector* makeOneStep(const MyVector& E, const MyVector& B);
	MyVector* getP();
	~TParticle();
};

