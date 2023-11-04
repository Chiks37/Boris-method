#include "TParticle.h"

TParticle::TParticle(double _x, double _y, double _z,
	double _px, double _py, double _pz,
	/*double _gamma, */double _q/*, double _velocity*/) :
	/*gamma(_gamma), */q(_q)//, v_new(_velocity)
{
	r_old.push_back(_x);
	r_old.push_back(_y);
	r_old.push_back(_z);
	p_old.push_back(_px);
	p_old.push_back(_py);
	p_old.push_back(_pz);
	r_new.resize(3);
	v_new.resize(3);
	t.resize(3);
	pMinus.resize(3);
	pPlus.resize(3);
	p_new.resize(3);
	s.resize(3);
	pDeriv.resize(3);
}

double TParticle::absValue(std::vector<double> vec)
{
	double sqrSum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	return pow(sqrSum, 1/2);
}

std::vector<double> TParticle::vecMul(std::vector<double> vec1, std::vector<double> vec2)
{
	std::vector<double> res;
	res.resize(3);
	double a = vec1[2] * vec2[1];
	double b = vec1[1] * vec2[2];
	res[0] =  b - a;
	res[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
	return res;
}

void TParticle::makeNewPMinus(std::vector<double> E)
{
	for (int i = 0; i < 3; i++) {
		pMinus[i] = p_old[i] + q * E[i] * (delta_t / 2);
	}
}

void TParticle::makeNewGammaOld()
{
	gamma_old = sqrt(1 + pow(absValue(p_old) / (m * c), 2));
}

void TParticle::makeNewt(std::vector<double> B)
{
	for (int i = 0; i < 3; i++) {
		t[i] = (q * B[i] * delta_t) / (gamma_old * m * c * 2);
	}
}

void TParticle::makeNewS()
{
	for (int i = 0; i < 3; i++) {
		s[i] = (2 * t[i]) / (1 + pow(absValue(t), 2));
	}
}

void TParticle::makeNewPDeriv()
{
	std::vector<double> mul = vecMul(pMinus, t);
	for (int i = 0; i < 3; i++) {
		pDeriv[i] = pMinus[i] + mul[i];
	}
}

void TParticle::makeNewPPlus()
{
	//std::vector<double> res;
	std::vector<double> mul = vecMul(pDeriv, s);
	for (int i = 0; i < 3; i++) {
		pPlus[i] = pMinus[i] + mul[i];
	}
}

void TParticle::makeNewPNew(std::vector<double> E)
{
	for (int i = 0; i < 3; i++) {
		p_new[i] = pPlus[i] + q * E[i] * delta_t / 2;
	}
}

void TParticle::makeNewGammaNew()
{
	gamma_new = sqrt(1 + pow(absValue(p_old) / (m * c), 2));
}

void TParticle::makeNewVNew()
{
	for (int i = 0; i < 3; i++) {
		v_new[i] = p_new[i] / (gamma_new * m);
	}
}

void TParticle::makeNewRNew()
{
	for (int i = 0; i < 3; i++) {
		r_new[i] = r_old[i] + v_new[i] * delta_t;
	}
}

void TParticle::updateAllAndPrint(std::vector<double> E, std::vector<double> B)
{
	makeNewPMinus(E);
	makeNewGammaOld();
	makeNewt(B);
	makeNewS();
	makeNewPDeriv();
	makeNewPPlus();
	makeNewPNew(E);
	makeNewGammaNew();
	makeNewVNew();
	makeNewRNew();
	std::cout << "Old particle coordinates:\n x: " << p_old[0] << "\n y: " << p_old[1] << "\n z: " << p_old[2] << "\n";
	std::cout << "New particle coordinates:\n x: " << p_new[0] << "\n y: " << p_new[1] << "\n z: " << p_new[2] << "\n";
}

