#include <iostream>
#include <vector>
#pragma once
class TParticle
{
	//Необходимость использование double или float требует дальнейшего исследования
	//позиция и момент - вектора трёх элементов!
	std::vector<double> position; // Particle position (x, y, z)
	std::vector<double> p_old; // Particle momentum (px, py, pz)
	double weight = 1; // Particle weight
	double gamma_old; // Particle γ-factor
	double q; // Particle type
	double velocity;
	double delta_t = 1; //1 second
	std::vector<double> pMinus;
	std::vector<double> pPLus;
	double absValue(std::vector<double> vec);
public:
	TParticle(double _x, double _y, double _z,
		double _px, double _py, double _pz,
		/*double _gamma, */double _q, double _velocity);
	void makeNewPMinus(std::vector<double> E);
	void makeNewGammaOld();
};

