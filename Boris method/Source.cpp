#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <omp.h>
#include <fstream>
#include "TParticle.h"
#include <math.h>
#include <cmath>

int NUM_THREADS = 1;

double rnd() {
	return rand() % 1000;
}

void calculate(TParticle* parts, const int& partsCount, const int& iterCount, const MyVector& E, const MyVector& B, double dt = 1e-12)
{
	double m = 9.10958215e-28;
	double q = -4.803e-10;
	const double mcbased = m * c;
	const double qdthbased = q * dt * 0.5;

	for (int j = 0; j < iterCount; j++)
	{
		//#pragma omp simd
		//#pragma ivdep
		//#pragma vector always
		//#pragma novector

		//#pragma parallel for

		omp_set_num_threads(NUM_THREADS);
		//#pragma omp parallel for
					//#pragma omp for simd
		//#pragma novector
#pragma omp parallel for simd 
//#pragma omp for nowait
		for (int i = 0; i < partsCount; i++)
		{
			double mc = mcbased;
			double qdth = qdthbased;
			double delta_t = dt;

			double temp_pMin = qdth;
			MyVector pMinus = parts[i].p_old + E * temp_pMin;

			double temp1_go = parts[i].p_old.absValue() / mc;
			double temp_go = temp1_go * temp1_go;
			double gamma_old = sqrt(1.0 + temp_go);

			MyVector t = (B * qdth) / (gamma_old * mc);

			double temp_s1 = t.absValue();
			double temp_s = temp_s1 * temp_s1;
			MyVector s = t * 2.0 / (1.0 + temp_s);

			MyVector pDeriv = pMinus + pMinus.vecMul(t);
			MyVector pPlus = pMinus + pDeriv.vecMul(s);
			MyVector p_new = pPlus + E * qdth;

			double temp1_gn = p_new.absValue();
			double temp2_gn = mc;
			double temp_gn = temp1_gn / temp2_gn;
			double tempSqr_gn = temp_gn * temp_gn;
			double gamma_new = sqrt(1.0 + tempSqr_gn);

			MyVector v_new = p_new / (gamma_new * m);
			MyVector r_new = parts[i].r_old + v_new * delta_t;

			parts[i].p_old = p_new;
			parts[i].r_old = r_new;
		}
	}
}

void makeOneStep(TParticle& part, const int& iterCount, const MyVector& E, const MyVector& B, double dt)
{
	TParticle* dynamicPart = new TParticle(part.r_old[0], part.r_old[1], part.r_old[2], part.p_old[0], part.p_old[1], part.p_old[2]);
	calculate(dynamicPart, 1, iterCount, E, B, dt);
	part.p_old = dynamicPart->p_old;
	part.r_old = dynamicPart->r_old;
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
	double delta_t = dt;

	makeOneStep(part, stepsCount, E, B, dt);

	MyVector rResult = part.r_old;
	MyVector pResult = part.p_old;

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
	const double dt = M_PI * me * c * sqrt(1 + pow(px / (me * c), 2)) / (qe * (-1) * Bz * stepsCount);

	MyVector E;
	E[0] = 0;
	E[1] = 0;
	E[2] = 0;

	MyVector B;
	B[0] = 0;
	B[1] = 0;
	B[2] = Bz;

	TParticle part(0, 0, 0, px, 0, 0);
	double delta_t = dt;

	makeOneStep(part, stepsCount, E, B, dt);

	MyVector rResult = part.r_old;
	MyVector pResult = part.p_old;

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

		if (tmpStream >> partsCount)
		{
			if (argc >= 4)
			{
				std::stringstream tmpStream(argv[3]);
				int tmp;
				if (tmpStream >> tmp && tmp > 0 && tmp <= 32)
				{
					NUM_THREADS = tmp;
				}
				else
				{
					std::cout << "Wrong threads num. Default = 1.\n";
				}
			}
			else
			{
				std::cout << "Default threads num = 1. \n";
			}
			std::cout << "\n\nNUM THREADS:" << NUM_THREADS << '\n';

			TParticle* parts = new TParticle[partsCount];
#pragma omp parallel for
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
		else
		{
			std::cout << "Wrong argument.\n";
		}
	}
	else
	{
		std::cout << "Need to pass argument.\n";
	}
	return 0;
}