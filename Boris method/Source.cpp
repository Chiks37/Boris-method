#include <iostream>
#include <sstream>
#include <omp.h>
#include <fstream>
#include "TParticle.h"

double rnd() {
	return rand() % 1000;
}

void main(int argc, char* argv[]) { // ¬сего 1 аргумент - количество временных шагов
	if (argc >= 2) {

		std::vector<double> tmp = { 0, 0, 0 };
		MyVector E;
		E = tmp;

		tmp = { 0.1, 0, 0 };
		MyVector B;
		B = tmp;

		std::stringstream tmpStream(argv[1]);
		double stepsCount;

		if (tmpStream >> stepsCount) {

			TParticle part(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());

			std::ofstream fout("partData", std::ios_base::app);
			fout.clear();
			MyVector r;
			double time = 0;
			for (int i = 0; i < stepsCount; i++) {

				double start = omp_get_wtime();
				r = part.makeOneStep(E, B);
				double end = omp_get_wtime();
				time += end - start;



				std::cout << r[0] << ' ' << r[1] << ' ' << r[2] << std::endl;
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