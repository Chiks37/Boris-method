#include "TParticle.h"

TParticle::TParticle() : q(1), gamma_old(0), gamma_new(0)
{
	r_old[0] = 0;
	r_old[1] = 0;
	r_old[2] = 0;
	p_old[0] = 0;
	p_old[1] = 0;
	p_old[2] = 0;
}

TParticle::TParticle(double _x, double _y, double _z,
	double _px, double _py, double _pz, double _q) :	q(_q), gamma_old(0), gamma_new(0)
{
	r_old[0] = _x;
	r_old[1] = _y;
	r_old[2] = _z;
	p_old[0] = _px;
	p_old[1] = _py;
	p_old[2] = _pz;
}

void TParticle::updateAllAndPrint(const MyVector& E, const MyVector& B, const size_t& K)
{
	pMinus = p_old + E * q * delta_t / 2;
	gamma_old = sqrt(1 + pow(p_old.absValue() / (m * c), 2));
	t = (B * q * delta_t) / (gamma_old * m * c * 2) * K;
	s = t * 2 / (1 + pow(t.absValue(), 2));
	pDeriv = pMinus + pMinus.vecMul(t);
	pPlus = pMinus + pDeriv.vecMul(s);
	p_new = pPlus + E * q * delta_t / 2;
	gamma_new = sqrt(1 + pow(p_old.absValue() / (m * c), 2));
	v_new = p_new / (gamma_new * m);
	r_new = r_old + v_new * delta_t;
	std::ofstream fout("partData", std::ios_base::app);
	if (K == 1) {
		fout << r_old[0] << ' ' << r_old[1] << ' ' << r_old[2] << std::endl;
	}
	fout /*<< r_old[0] << ' ' << r_old[1] << ' ' << r_old[2] << ' ' */<< r_new[0] << ' ' << r_new[1] << ' ' << r_new[2] << std::endl;
	fout.close();
	//std::cout << "Old particle coordinates:\n x: " << r_old[0] << "\n y: " << r_old[1] << "\n z: " << r_old[2] << "\n";
	//std::cout << "New particle coordinates:\n x: " << r_new[0] << "\n y: " << r_new[1] << "\n z: " << r_new[2] << "\n";
	//std::cout << "Old particle coordinates:\n x: " << p_old[0] << "\n y: " << p_old[1] << "\n z: " << p_old[2] << "\n";
	//std::cout << "New particle coordinates:\n x: " << p_new[0] << "\n y: " << p_new[1] << "\n z: " << p_new[2] << "\n";
}

