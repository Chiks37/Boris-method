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

bool relAccelInStatFieldTest() {

	const double qe = -4.803e-10;
	const double me = 9.10958215e-28;
	const double Ex = 2;
	const int stepsCount = 100;
	const double dt = me * c / (qe * Ex * stepsCount * 1e6);
	//const double v0 = 0;

	std::vector<double> tmp = { Ex, 0, 0 };
	MyVector E;
	E = tmp;

	tmp = { 0, 0, 0 };
	MyVector B;
	B = tmp;

	TParticle part(0, 0, 0, 0, 0, 0);
	part.m = me;
	part.q = qe;
	part.delta_t = dt;

	for (int i = 0; i < stepsCount - 1; i++) {
		part.makeOneStep(E, B);
	}
	MyVector rResult = part.makeOneStep(E, B);
	MyVector pResult = part.getP();

	tmp = { me * c * c / (2 * qe * Ex * 1e12), 0, 0 };
	MyVector rAnalitic = tmp;

	tmp = { me * c / 1e6, 0, 0 };
	MyVector pAnalitic = tmp;

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

	std::vector<double> tmp = { 0, 0, 0 };
	MyVector E;
	E = tmp;

	tmp = { 0, 0, Bz };
	MyVector B;
	B = tmp;

	TParticle part(0, 0, 0, px, 0, 0);
	part.m = me;
	part.q = qe;
	part.delta_t = dt;

	for (int i = 0; i < stepsCount - 1; i++) {
		part.makeOneStep(E, B);
	}
	MyVector rResult = part.makeOneStep(E, B);
	MyVector pResult = part.getP();

	tmp = { 0, -2 * px * c / (qe * Bz), 0 };
	MyVector rAnalitic = tmp;

	tmp = { -px, 0, 0 };
	MyVector pAnalitic = tmp;

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

void main(int argc, char* argv[]) { // 2 ��������� - ���������� ��������� � ��������
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

		std::vector<double> tmp = { 2, 2, 2 };
		MyVector E(tmp);

		tmp = { 0.1, 0.1, 0.1 };
		MyVector B(tmp);

		std::stringstream tmpStream(argv[1]);
		double partsCount;

		if (tmpStream >> partsCount) {
			TParticle* parts = new TParticle[partsCount];
			for (int i = 0; i < partsCount; i++) {
				parts[i] = TParticle(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
			}

			double start = omp_get_wtime();
			for (int j = 0; j < iterCount; j++) {
				for (int i = 0; i < partsCount; i++) {
					parts[i].makeOneStep(E, B);
				}
			}
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
}