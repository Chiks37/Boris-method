#include "TParticle.h"

TParticle::TParticle(double _x, double _y, double _z,
	double _px, double _py, double _pz,
	/*double _gamma, */double _q, double _velocity) : 
	/*gamma(_gamma), */q(_q), velocity(_velocity)
{
	position.push_back(_x);
	position.push_back(_y);
	position.push_back(_z);
	p_old.push_back(_px);
	p_old.push_back(_py);
	p_old.push_back(_pz);
}

double TParticle::absValue(std::vector<double> vec)
{
	double sqrSum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	return pow(sqrSum, 1/2);
}

void TParticle::makeNewPMinus(std::vector<double> E)
{
	for (int i = 0; i < 3; i++) {
		pMinus[i] = p_old[i] + q * E[i] * (delta_t / 2);
	}
}

void TParticle::makeNewGammaOld()
{

}
