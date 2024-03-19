#include "TParticle.h"

TParticle::TParticle() //: gamma_old(0), gamma_new(0)
{
	r_old[0] = 0;
	r_old[1] = 0;
	r_old[2] = 0;

	p_old[0] = 0;
	p_old[1] = 0;
	p_old[2] = 0;
}

TParticle::TParticle(double _x, double _y, double _z,
	double _px, double _py, double _pz) //: gamma_old(0), gamma_new(0)
{
	r_old[0] = _x;
	r_old[1] = _y;
	r_old[2] = _z;

	p_old[0] = _px;
	p_old[1] = _py;
	p_old[2] = _pz;
}

//MyVector TParticle::makeOneStep(const MyVector& E, const MyVector& B)
//{
//	pMinus = p_old + E * q * delta_t / 2;
//	gamma_old = sqrt(1 + (p_old.absValue() / (m * c)) * (p_old.absValue() / (m * c)));
//	t = (B * q * delta_t) / (gamma_old * m * c * 2);
//	s = t * 2 / (1 + (t.absValue()) * (t.absValue()));
//	pDeriv = pMinus + pMinus.vecMul(t);
//	pPlus = pMinus + pDeriv.vecMul(s);
//	p_new = pPlus + E * q * delta_t / 2;
//	gamma_new = sqrt(1 + (p_new.absValue() / (m * c)) * (p_new.absValue() / (m * c)));
//	v_new = p_new / (gamma_new * m);
//	r_new = r_old + v_new * delta_t;
//
//	p_old = p_new;
//	r_old = r_new;
//	v_old = v_new;
//
//	return r_new;
//}

MyVector TParticle::getP()
{
	return p_old;
}

