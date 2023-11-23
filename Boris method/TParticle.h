#include <iostream>
#include <vector>
#include "MyVector.h"
#include <fstream>
#pragma once
#define c 299792458
class TParticle
{
	//Необходимость использование double или float требует дальнейшего исследования
	double m = 1; // Particle weight
	double gamma_old; // Particle γ-factor
	double gamma_new;
	double q; // Particle type
	double delta_t = 1000000000; //1 second
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
		double _px, double _py, double _pz, double _q);

	void updateAllAndPrint(const MyVector& E, const MyVector& B, const size_t& K);
};

