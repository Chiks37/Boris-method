#include <iostream>
#include <vector>
#pragma once
#define c 299792458
class TParticle
{
	//Необходимость использование double или float требует дальнейшего исследования
	double m = 1; // Particle weight
	double gamma_old; // Particle γ-factor
	double gamma_new;
	double q; // Particle type
	double delta_t = 1; //1 second
	std::vector<double> r_old; // Particle position (x, y, z)
	std::vector<double> r_new;
	std::vector<double> v_new;
	std::vector<double> t;
	std::vector<double> p_old; // Particle momentum (px, py, pz)
	std::vector<double> pMinus;
	std::vector<double> pPlus;
	std::vector<double> p_new;
	std::vector<double> s;
	std::vector<double> pDeriv;
	double absValue(std::vector<double> vec);
	std::vector<double> vecMul(std::vector<double> vec1, std::vector<double> vec2);
public:

	TParticle(double _x, double _y, double _z,
		double _px, double _py, double _pz,
		/*double _gamma, */double _q/*, double _velocity*/);

	void makeNewPMinus(std::vector<double> E);
	void makeNewGammaOld();
	void makeNewt(std::vector<double> B);
	void makeNewS();
	void makeNewPDeriv();
	void makeNewPPlus();
	void makeNewPNew(std::vector<double> E);
	void makeNewGammaNew();
	void makeNewVNew();
	void makeNewRNew();
	void updateAllAndPrint(std::vector<double> E, std::vector<double> B);
};

