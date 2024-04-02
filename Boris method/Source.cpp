#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <fstream>
#include "TParticle.h"
#include <math.h>
#include <cmath>
#include <chrono>
#include <riscv-vector.h>

double rnd() {
	return rand() % 1000;
}

void calculate(TParticle* parts, const int& partsCount, const int& iterCount, const MyVector& E, const MyVector& B, double dt = 1e-12)
{
	double m = 9.10958215e-28;
	double q = -4.803e-10;

	double mcbased = m * c;
	double qdtbased = q * dt * 0.5;

	int vl = vsetvlmax_e32m4();
	vfloat32m4_t mcbasedv = vfmv_v_f_f32m4( mcbased, vl);

	vfloat32m4_t qdthbasedv = vfmv_v_f_f32m4(qdtbased, vl);

	vfloat32m4_t dtv = vfmv_v_f_f32m4(dt, vl);

	for (int j = 0; j < iterCount; j++)
	{
		for (int i = 0; i < partsCount; i++)
		{
			vfloat32m4_t mcv = mcbasedv;
			vfloat32m4_t qdthv = qdthbasedv;
			vfloat32m4_t delta_tv = dtv;

			//MyVector pMinus = parts[i].p_old + E * qdthv;
			vfloat32m4_t p_old1v = parts[i].p_old[0];
			vfloat32m4_t pMinus1v = /////////////////////// EDIT FROM HERE

			double temp1_go = parts[i].p_old.absValue() / mcv;
			double temp_go = temp1_go * temp1_go;
			double gamma_old = sqrt(1.0 + temp_go);

			MyVector t = (B * qdthv) / (gamma_old * mcv);

			double temp_s1 = t.absValue();
			double temp_s = temp_s1 * temp_s1;
			MyVector s = t * 2.0 / (1.0 + temp_s);

			MyVector pDeriv = pMinus + pMinus.vecMul(t);
			MyVector pPlus = pMinus + pDeriv.vecMul(s);
			MyVector p_new = pPlus + E * qdthv;

			double temp1_gn = p_new.absValue();
			double temp2_gn = mcv;
			double temp_gn = temp1_gn / temp2_gn;
			double tempSqr_gn = temp_gn * temp_gn;
			double gamma_new = sqrt(1.0 + tempSqr_gn);

			MyVector v_new = p_new / (gamma_new * m);
			MyVector r_new = parts[i].r_old + v_new * delta_tv;

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

int main(int argc, char* argv[]) { // 2 аргумента - количество элементов и итераций
	relAccelInStatFieldTest();
	osciInStaticMagFieldTest();
	if (argc >= 2) 
	{
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
			TParticle* parts = new TParticle[partsCount];
			for (int i = 0; i < partsCount; i++) {
				parts[i] = TParticle(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
			}

			auto start = std::chrono::high_resolution_clock::now(); // фиксация времени начала выполнения функции
			calculate(parts, partsCount, iterCount, E, B);
			auto end = std::chrono::high_resolution_clock::now(); // фиксация времени окончания выполнения функции

			std::chrono::duration<double> time = end - start; // вычисление продолжительности выполнения функции
			std::cout << "Execution time: " << time.count() << " seconds.\n";

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