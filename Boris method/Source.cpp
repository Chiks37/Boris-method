#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <omp.h>
#include <fstream>
#include "TParticle.h"
#include <math.h>
#include <risc-v.h>

//double rnd() {
//	return rand() % 1000;
//}

void calculate(TParticle& parts, int partsCount, int iterCount, const MyVector& E, const MyVector& B, double dt = 1e-12)
{
	size_t vl = vsetvl_e32m1(parts.partsCount);

	double m = 9.10958215e-28;  // use more exact constants
	vfloat32m1_t mV = vfmv_v_f_f32m1(m, vl);

	double* gamma_old = new double[parts.partsCount];
	vfloat32m1_t goV = vle_v_f32m1(gamma_old, vl);

	double* gamma_new = new double[parts.partsCount];
	vfloat32m1_t gnV = vle_v_f32m1(gamma_old, vl);

	double q = -4.803e-10;
	vfloat32m1_t qV = vfmv_v_f_f32m1(q, vl);

	double delta_t = dt;
	vfloat32m1_t dtV = vfmv_v_f_f32m1(dt, vl);


	//MyVector* r_new = new MyVector[parts.partsCount];
	double* r_new0 = new double[parts.partsCount];
	vfloat32m1_t rn0V = vle_v_f32m1(r_new0, vl);

	double* r_new1 = new double[parts.partsCount];
	vfloat32m1_t rn1V = vle_v_f32m1(r_new1, vl);

	double* r_new2 = new double[parts.partsCount];
	vfloat32m1_t rn2V = vle_v_f32m1(r_new2, vl);


	//MyVector* v_old = new MyVector[parts.partsCount];
	double* v_old0 = new double[parts.partsCount];
	vfloat32m1_t vo0V = vle_v_f32m1(v_old0, vl);

	double* v_old1 = new double[parts.partsCount];
	vfloat32m1_t vo1V = vle_v_f32m1(v_old1, vl);

	double* v_old2 = new double[parts.partsCount];
	vfloat32m1_t vo2V = vle_v_f32m1(v_old2, vl);


	//MyVector* v_new = new MyVector[parts.partsCount];
	double* v_new0 = new double[parts.partsCount];
	vfloat32m1_t vn0V = vle_v_f32m1(v_new0, vl);

	double* v_new1 = new double[parts.partsCount];
	vfloat32m1_t vn1V = vle_v_f32m1(v_new1, vl);

	double* v_new2 = new double[parts.partsCount];
	vfloat32m1_t vn2V = vle_v_f32m1(v_new2, vl);

	//MyVector* t = new MyVector[parts.partsCount];
	double* t0 = new double[parts.partsCount];
	vfloat32m1_t t0V = vle_v_f32m1(t0, vl);

	double* t1 = new double[parts.partsCount];
	vfloat32m1_t t1V = vle_v_f32m1(t1, vl);

	double* t2 = new double[parts.partsCount];
	vfloat32m1_t t2V = vle_v_f32m1(t2, vl);


	//MyVector* pMinus = new MyVector[parts.partsCount];
	double* pm0 = new double[parts.partsCount];
	vfloat32m1_t pm0V = vle_v_f32m1(pm0, vl);

	double* pm1 = new double[parts.partsCount];
	vfloat32m1_t pm1V = vle_v_f32m1(pm1, vl);

	double* pm2 = new double[parts.partsCount];
	vfloat32m1_t pm2V = vle_v_f32m1(pm2, vl);


	//MyVector* pPlus = new MyVector[parts.partsCount];
	double* pp0 = new double[parts.partsCount];
	vfloat32m1_t pp0V = vle_v_f32m1(pp0, vl);

	double* pp1 = new double[parts.partsCount];
	vfloat32m1_t pp1V = vle_v_f32m1(pp1, vl);

	double* pp0 = new double[parts.partsCount];
	vfloat32m1_t pp2V = vle_v_f32m1(pp2, vl);


	//MyVector* p_new = new MyVector[parts.partsCount];
	double* pn0 = new double[parts.partsCount];
	vfloat32m1_t pn0V = vle_v_f32m1(pn0, vl);

	double* pn1 = new double[parts.partsCount];
	vfloat32m1_t pn1V = vle_v_f32m1(pn1, vl);

	double* pn2 = new double[parts.partsCount];
	vfloat32m1_t pn2V = vle_v_f32m1(pn2, vl);


	MyVector* s = new MyVector[parts.partsCount];
	MyVector* pDeriv = new MyVector[parts.partsCount];

	double mcbased = m * c;
	double qdthbased = q * delta_t * 0.5;

	for (int j = 0; j < iterCount; j++)
	{
		for (int i = 0; i < partsCount; i++)
		{
			double mc = mcbased;
			double qdth = qdthbased;

			pMinus[i] = parts.p_old[i] + E * qdth;
			gamma_old[i] = sqrt(1.0 + (parts.p_old[i].absValue() / mc) * (parts.p_old[i].absValue() / mc));
			t[i] = (B * qdth) / (gamma_old[i] * mc);
			s[i] = t[i] * 2.0 / (1.0 + (t[i].absValue()) * (t[i].absValue()));
			pDeriv[i] = pMinus[i] + pMinus[i].vecMul(t[i]);
			pPlus[i] = pMinus[i] + pDeriv[i].vecMul(s[i]);
			p_new[i] = pPlus[i] + E * qdth;
			gamma_new[i] = sqrt(1 + (p_new[i].absValue() / mc) * (p_new[i].absValue() / mc));
			v_new[i] = p_new[i] / (gamma_new[i] * m);
			r_new[i] = parts.r_old[i] + v_new[i] * delta_t;

			parts.p_old[i] = p_new[i];
			parts.r_old[i] = r_new[i];

		/*	if (i == 0 && j == 0)
			{
				std::cout << pMinus[i] << '\n' << gamma_old[i] << '\n' << t[i] << '\n' << s[i] << '\n' << pDeriv[i] << '\n' << pPlus[i] << '\n' << p_new[i] << '\n' << gamma_new[i] << '\n' << v_new[i] << '\n' << r_new[i] << '\n' << parts.p_old[i] << '\n' << parts.r_old[i] << '\n';
			}*/
		}
		//parts.makeOneStep(E, B);
	}
	delete[] gamma_old;
	delete[] gamma_new;
	delete[] r_new;
	delete[] v_old;
	delete[] v_new;
	delete[] t;
	delete[] pMinus;
	delete[] pPlus;
	delete[] p_new;
	delete[] s;
	delete[] pDeriv;
}

void makeOneStep(TParticle& part, const int& iterCount, const MyVector& E, const MyVector& B, double dt)
{
	calculate(part, 1, iterCount, E, B, dt);
}

bool relAccelInStatFieldTest() {

	const double qe = -4.803e-10;
	const double me = 9.10958215e-28;
	const double Ex = 2;
	const int stepsCount = 100;
	const double dt = me * c / (qe * Ex * stepsCount * 1e6);

	MyVector E;
	E[0] = Ex;
	E[1] = 0;
	E[2] = 0;

	MyVector B;
	B[0] = 0;
	B[1] = 0;
	B[2] = 0;

	TParticle part(0, 0, 0, 0, 0, 0, 1);
	double m = me;
	double q = qe;
	//double delta_t = dt;
	double qdthbased = qe * dt * 0.5;
	double mcbased = me * c;

	makeOneStep(part, stepsCount, E, B, dt);

	MyVector rResult = part.r_old[0];
	MyVector pResult = part.getP()[0];

	MyVector rAnalitic;
	rAnalitic[0] = me * c * c / (2 * qe * Ex * 1e12);
	rAnalitic[1] = 0;
	rAnalitic[2] = 0;

	MyVector pAnalitic;
	pAnalitic[0] = me * c / 1e6;
	pAnalitic[1] = 0;
	pAnalitic[2] = 0;

	std::cout << "Relativistic acceleration in statick field test:\nCalculated coordinates: ";
	std::cout << rResult[0] << ' ' << rResult[1] << ' ' << rResult[2] << std::endl;
	std::cout << "Calculated impulse:     ";
	std::cout << pResult[0] << ' ' << pResult[1] << ' ' << pResult[2] << std::endl;
	std::cout << "Analitic coordinates:   ";
	std::cout << rAnalitic[0] << ' ' << rAnalitic[1] << ' ' << rAnalitic[2] << std::endl;
	std::cout << "Analitic impulse:       ";
	std::cout << pAnalitic[0] << ' ' << pAnalitic[1] << ' ' << pAnalitic[2] << "\n\n";
	return true;
}

bool osciInStaticMagFieldTest() {

	const double qe = -4.803e-10;
	const double me = 9.10958215e-28;
	const double Bz = 0.1;
	const int stepsCount = 100;
	const double v0 = 3e8;
	const double px = me * v0; //recheck
	const double dt = M_PI * me * c * sqrt(1 + pow(px / (me * c), 2)) / (abs(qe) * Bz * stepsCount);

	MyVector E;
	E[0] = 0;
	E[1] = 0;
	E[2] = 0;

	MyVector B;
	B[0] = 0;
	B[1] = 0;
	B[2] = Bz;

	TParticle part(0, 0, 0, px, 0, 0, 1);
	double m = me;
	double q = qe;
	//double delta_t = dt;
	double qdthbased = qe * dt * 0.5;
	double mcbased = me * c;

	makeOneStep(part, stepsCount, E, B, dt);

	MyVector rResult = part.r_old[0];
	MyVector pResult = part.getP()[0];

	MyVector rAnalitic;
	rAnalitic[0] = 0;
	rAnalitic[1] = -2 * px * c / (qe * Bz);
	rAnalitic[2] = 0;

	MyVector pAnalitic;
	pAnalitic[0] = -px;
	pAnalitic[1] = 0;
	pAnalitic[2] = 0;

	std::cout << "Oscilation in a static magnet field test:\nCalculated coordinates: ";
	std::cout << rResult[0] << ' ' << rResult[1] << ' ' << rResult[2] << std::endl;
	std::cout << "Calculated impulse:     ";
	std::cout << pResult[0] << ' ' << pResult[1] << ' ' << pResult[2] << std::endl;
	std::cout << "Analitic coordinates:   ";
	std::cout << rAnalitic[0] << ' ' << rAnalitic[1] << ' ' << rAnalitic[2] << std::endl;
	std::cout << "Analitic impulse:       ";
	std::cout << pAnalitic[0] << ' ' << pAnalitic[1] << ' ' << pAnalitic[2] << "\n\n";
	return true;
}

// TEST CONFIGURATION:
// non-relativistic case, v << c, gamma->1
// Ex = Ey = Bx = By = 0, Ez=E, Bz=B
// v(t) = (vx(t), vy(t), vz(t)), v(0) = (0, v0, 0)
// r(t) = (rx(t), ry(t), rz(t)), r(0) = (-v0/omega, 0, 0)
// omega = qB/(mc)
//
// ANALYTICAL SOLUTION:
// v(t) = (v0*sin(omega*t), v0*cos(omega*t), qE/m*t)
// r(t) = (-v0/omega*cos(omega*t), v0/omega*sin(omega*t), qE/m*t*t/2)
// x, y axes -> motion in a circle
// z axis -> uniformly accelerated motion

int main(int argc, char* argv[]) { // 2 аргумента - количество элементов и итераций
	relAccelInStatFieldTest();
	osciInStaticMagFieldTest();
	if (argc >= 2) {

		//const double qe = -4.803e-10;
		//const double me = 9.10958215e-28;
		//const double Bz = 0.1;
		//const double Ez = 0.00001;  // small E to avoid relativistic case
		//const double omega = qe * Bz / me / c;
		//const double v0 = omega;  // v0 << c, non-relativistic case
		//const double gamma0 = 1.0 / sqrt(1.0 - v0 * v0 / (c * c));
		//const double dt = 1e-13*omega;
		//const double T = 2 * 3.1415 / omega;  // period, dt should be less than T

		int iterCount = 100;
		if (argc >= 3) {
			std::stringstream tmpStream(argv[2]);
			if (!(tmpStream >> iterCount)) {
				std::cout << "Wrong second argument.\n";
			}
		}

		MyVector E;
		E[0] = 2;
		E[1] = 2;
		E[2] = 2;

		MyVector B;
		B[0] = 0.1;
		B[1] = 0.1;
		B[2] = 0.1;

		std::stringstream tmpStream(argv[1]);
		int  partsCount;

		if (tmpStream >> partsCount && partsCount != 0) {
			TParticle parts(partsCount);

			double start = omp_get_wtime();
			calculate(parts, partsCount, iterCount, E, B);
			double end = omp_get_wtime();

			double time = end - start;
			std::cout << "Execution time: " << time << " seconds.\n";
		}
		else {
			std::cout << "Wrong argument.\n";
		}
	}
	else {
		std::cout << "Need to pass argument.\n";
	}
	return 0;
}