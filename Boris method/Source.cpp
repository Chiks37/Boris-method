#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <omp.h>
#include <fstream>
#include "TParticle.h"
#include <math.h>

double rnd() {
	return rand() % 1000;
}

void calculate(TParticle* parts, const int partsCount, const int iterCount, const MyVector& E, const MyVector& B)
{
	for (int j = 0; j < iterCount; j++)
	{
#pragma omp simd
//#pragma ivdep
//#pragma vector always
//#pragma novector
		for (int i = 0; i < partsCount; i++)
		{
			parts[i].pMinus = parts[i].p_old + E * parts[i].q * parts[i].delta_t * 0.5;
			parts[i].gamma_old = sqrt(1.0 + (parts[i].p_old.absValue() / (parts[i].m * c)) * (parts[i].p_old.absValue() / (parts[i].m * c)));
			parts[i].t = (B * parts[i].q * parts[i].delta_t) / (parts[i].gamma_old * parts[i].m * c * 2.0);
			parts[i].s = parts[i].t * 2.0 / (1.0 + (parts[i].t.absValue()) * (parts[i].t.absValue()));
			parts[i].pDeriv = parts[i].pMinus + parts[i].pMinus.vecMul(parts[i].t);
			parts[i].pPlus = parts[i].pMinus + parts[i].pDeriv.vecMul(parts[i].s);
			parts[i].p_new = parts[i].pPlus + E * parts[i].q * parts[i].delta_t * 0.5;
			parts[i].gamma_new = sqrt(1.0 + (parts[i].p_new.absValue() / (parts[i].m * c)) * (parts[i].p_new.absValue() / (parts[i].m * c)));
			parts[i].v_new = parts[i].p_new / (parts[i].gamma_new * parts[i].m);
			parts[i].r_new = parts[i].r_old + parts[i].v_new * parts[i].delta_t;

			parts[i].p_old = parts[i].p_new;
			parts[i].r_old = parts[i].r_new;
			parts[i].v_old = parts[i].v_new;
			//parts[i].makeOneStep(E, B);
		}
	}
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

	TParticle part(0, 0, 0, 0, 0, 0);
	part.m = me;
	part.q = qe;
	part.delta_t = dt;

	for (int i = 0; i < stepsCount - 1; i++) {
		part.makeOneStep(E, B);
	}
	MyVector rResult = part.makeOneStep(E, B);
	MyVector pResult = part.getP();

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

	TParticle part(0, 0, 0, px, 0, 0);
	part.m = me;
	part.q = qe;
	part.delta_t = dt;

	for (int i = 0; i < stepsCount - 1; i++) {
		part.makeOneStep(E, B);
	}
	MyVector rResult = part.makeOneStep(E, B);
	MyVector pResult = part.getP();

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

		if (tmpStream >> partsCount) {
			TParticle* parts = new TParticle[partsCount];
			for (int i = 0; i < partsCount; i++) {
				parts[i] = TParticle(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
			}

			double start = omp_get_wtime();
			calculate(parts, partsCount, iterCount, E, B);
			double end = omp_get_wtime();

			double time = end - start;
			std::cout << "Execution time: " << time << " seconds.\n";

			delete[] parts;
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