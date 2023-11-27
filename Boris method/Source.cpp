#include <iostream>
#include <sstream>
#include <omp.h>
#include <fstream>
#include "TParticle.h"

double rnd() {
	return rand() % 1000;
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

void main(int argc, char* argv[]) { // ¬сего 1 аргумент - количество временных шагов
	if (argc >= 2) {

		const double qe = -4.803e-10;
		const double me = 9.10958215e-28;
		const double Bz = 0.1;
		const double Ez = 0.00001;  // small E to avoid relativistic case
		const double omega = qe * Bz / me / c;
		const double v0 = omega;  // v0 << c, non-relativistic case
		const double gamma0 = 1.0 / sqrt(1.0 - v0 * v0 / (c * c));
		const double dt = 1e-13*omega;
		const double T = 2 * 3.1415 / omega;  // period, dt should be less than T

		std::vector<double> tmp = { 0, 0, Ez };
		MyVector E;
		E = tmp;

		tmp = { 0, 0, Bz };
		MyVector B;
		B = tmp;

		std::stringstream tmpStream(argv[1]);
		double stepsCount;

		if (tmpStream >> stepsCount) {

			//TParticle part(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
			TParticle part(-v0 / omega, 0, 0, 0, gamma0*me*v0, 0);
			part.delta_t = dt;

			std::ofstream fout("partData");  // clear and open for writing

			MyVector r;
			double time = 0;
			for (int i = 0; i < stepsCount; i++) {

				double start = omp_get_wtime();
				r = part.makeOneStep(E, B);
				double end = omp_get_wtime();
				time += end - start;

				fout << r[0] << ' ' << r[1] << ' ' << r[2] << std::endl;
				
				//double t = (i + 1) * dt;
				//fout << -v0 / omega * cos(omega * t) << ' ' << v0 / omega * sin(omega * t) << ' ' << qe*Ez / me * t * t * 0.5 << std::endl;
			}
			fout.close();

			std::cout << "Execution time: " << time << " seconds.\n";

		}
		else {
			std::cout << "Wrong argument.\n";
		}
	}
	else {
		std::cout << "Need to pass argument.\n";
	}
	system("python pointsOutput.py");
}