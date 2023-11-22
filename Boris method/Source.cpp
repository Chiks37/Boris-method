#include <iostream>
#include <sstream>
#include "TParticle.h"
#include <omp.h>

double rnd() {
	return rand() % 1000;
}
void main(int argc, char* argv[]) {
	//srand(time(0));
	if (argc >= 2) {
		double timeSum = 0;

		int itCount = 1; //iterations count
		if (argc >= 3) {
			std::stringstream tmpStream(argv[2]);
			if (!(tmpStream >> itCount)) {
				itCount = 1;
				std::cout << "Wrong second argument.\n";
			}
		}

		for (int j = 0; j < itCount; j++) {
			std::vector<double> tmp = { 2, 2, 2 };
			MyVector E;
			E = tmp;

			tmp = { 0.1, 0.1, 0.1 };
			MyVector B;
			B = tmp;

			std::stringstream tmpStream(argv[1]);
			double count;
			if (tmpStream >> count) {
				TParticle* parts = new TParticle[count];
				for (int i = 0; i < count; i++) {
					parts[i] = TParticle(rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
				}
				//TParticle part(0, 0, 0, 1, 1, 1, 1);

				double start = omp_get_wtime();
				for (int i = 0; i < count; i++) {
					parts[i].updateAllAndPrint(E, B);
				}
				double end = omp_get_wtime();

				double time = end - start;
				//double time = ((double)end - start) / ((double)CLOCKS_PER_SEC);
				std::cout << "Execution time: " << time << " seconds.\n";
				timeSum += time;

				delete[] parts;
			}
			else {
				std::cout << "Wrong first argument.\n";
			}
		}
		std::cout << "  AVERAGE TIME: " << timeSum / itCount << " seconds.\n";
	}
	else {
		std::cout << "Need to pass argument.\n";
	}
	system("python pointsOutput.py");
}