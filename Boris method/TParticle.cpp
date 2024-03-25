#include "TParticle.h"

double rnd() {
	return rand() % 1000;
}

TParticle::TParticle(double _x, double _y, double _z,
	double _px, double _py, double _pz, int _partsCount) : partsCount(_partsCount)
{
	gamma_old = new double[partsCount];
	gamma_new = new double[partsCount];
	r_old = new MyVector[partsCount];
	r_new = new MyVector[partsCount];
	v_old = new MyVector[partsCount];
	v_new = new MyVector[partsCount];
	t = new MyVector[partsCount];
	p_old = new MyVector[partsCount];
	pMinus = new MyVector[partsCount];
	pPlus = new MyVector[partsCount];
	p_new = new MyVector[partsCount];
	s = new MyVector[partsCount];
	pDeriv = new MyVector[partsCount];

	for (int i = 0; i < partsCount; i++)
	{
		gamma_old[i] = 0;
		gamma_new[i] = 0;

		r_old[i][0] = _x;
		r_old[i][1] = _y;
		r_old[i][2] = _z;

		p_old[i][0] = _px;
		p_old[i][1] = _py;
		p_old[i][2] = _pz;
	}
}

TParticle::TParticle(int _partsCount) : partsCount(_partsCount)
{
	gamma_old = new double[partsCount];
	gamma_new = new double[partsCount];
	r_old = new MyVector[partsCount];
	r_new = new MyVector[partsCount];
	v_old = new MyVector[partsCount];
	v_new = new MyVector[partsCount];
	t = new MyVector[partsCount];
	p_old = new MyVector[partsCount];
	pMinus = new MyVector[partsCount];
	pPlus = new MyVector[partsCount];
	p_new = new MyVector[partsCount];
	s = new MyVector[partsCount];
	pDeriv = new MyVector[partsCount];

	for (int i = 0; i < partsCount; i++)
	{
		gamma_old[i] = 0;
		gamma_new[i] = 0;

		r_old[i][0] = rnd();
		r_old[i][1] = rnd();
		r_old[i][2] = rnd();

		p_old[i][0] = rnd();
		p_old[i][1] = rnd();
		p_old[i][2] = rnd();
	}
}

TParticle::TParticle(const TParticle& ob2) : partsCount(ob2.partsCount)
{
	gamma_old = new double[partsCount];
	gamma_new = new double[partsCount];
	r_old = new MyVector[partsCount];
	r_new = new MyVector[partsCount];
	v_old = new MyVector[partsCount];
	v_new = new MyVector[partsCount];
	t = new MyVector[partsCount];
	p_old = new MyVector[partsCount];
	pMinus = new MyVector[partsCount];
	pPlus = new MyVector[partsCount];
	p_new = new MyVector[partsCount];
	s = new MyVector[partsCount];
	pDeriv = new MyVector[partsCount];

	for (int i = 0; i < partsCount; i++)
	{
		gamma_old[i] = ob2.gamma_old[i];
		gamma_new[i] = ob2.gamma_new[i];

		r_old[i][0] = ob2.r_old[i][0];
		r_old[i][1] = ob2.r_old[i][1];
		r_old[i][2] = ob2.r_old[i][2];

		r_new[i][0] = ob2.r_new[i][0];
		r_new[i][1] = ob2.r_new[i][1];
		r_new[i][2] = ob2.r_new[i][2];

		v_old[i][0] = ob2.v_old[i][0];
		v_old[i][1] = ob2.v_old[i][1];
		v_old[i][2] = ob2.v_old[i][2];

		v_new[i][0] = ob2.v_new[i][0];
		v_new[i][1] = ob2.v_new[i][1];
		v_new[i][2] = ob2.v_new[i][2];

		t[i][0] = ob2.t[i][0];
		t[i][1] = ob2.t[i][1];
		t[i][2] = ob2.t[i][2];

		p_old[i][0] = ob2.p_old[i][0];
		p_old[i][1] = ob2.p_old[i][0];
		p_old[i][2] = ob2.p_old[i][0];

		p_new[i][0] = ob2.p_new[i][0];
		p_new[i][1] = ob2.p_new[i][1];
		p_new[i][2] = ob2.p_new[i][2];

		pMinus[i][0] = ob2.pMinus[i][0];
		pMinus[i][1] = ob2.pMinus[i][1];
		pMinus[i][2] = ob2.pMinus[i][2];

		pPlus[i][0] = ob2.pPlus[i][0];
		pPlus[i][1] = ob2.pPlus[i][1];
		pPlus[i][2] = ob2.pPlus[i][2];

		s[i][0] = ob2.s[i][0];
		s[i][1] = ob2.s[i][1];
		s[i][2] = ob2.s[i][2];

		pDeriv[i][0] = ob2.pDeriv[i][0];
		pDeriv[i][1] = ob2.pDeriv[i][1];
		pDeriv[i][2] = ob2.pDeriv[i][2];
	}

}

MyVector* TParticle::makeOneStep(const MyVector& E, const MyVector& B)
{
#pragma omp simd
//#pragma ivdep
//#pragma vector always
	for (int i = 0; i < partsCount; i++)
	{
		double mc = mcbased;
		double qdth = qdthbased;

		pMinus[i] = p_old[i] + E * qdth;
		gamma_old[i] = sqrt(1.0 + (p_old[i].absValue() / mc) * (p_old[i].absValue() / mc));
		t[i] = (B * qdth) / (gamma_old[i] * mc);
		s[i] = t[i] * 2.0 / (1.0 + (t[i].absValue()) * (t[i].absValue()));
		pDeriv[i] = pMinus[i] + pMinus[i].vecMul(t[i]);
		pPlus[i] = pMinus[i] + pDeriv[i].vecMul(s[i]);
		p_new[i] = pPlus[i] + E * qdth;
		gamma_new[i] = sqrt(1 + (p_new[i].absValue() / mc) * (p_new[i].absValue() / mc));
		v_new[i] = p_new[i] / (gamma_new[i] * m);
		r_new[i] = r_old[i] + v_new[i] * delta_t;

		p_old[i] = p_new[i];
		r_old[i] = r_new[i];
		v_old[i] = v_new[i];
	}

	return r_new;
}

MyVector* TParticle::getP()
{
	return p_old;
}

TParticle::~TParticle()
{
	delete[] gamma_old;
	delete[] gamma_new;
	delete[] r_old;
	delete[] r_new;
	delete[] v_old;
	delete[] v_new;
	delete[] t;
	delete[] p_old;
	delete[] pMinus;
	delete[] pPlus;
	delete[] p_new;
	delete[] s;
	delete[] pDeriv;
}

