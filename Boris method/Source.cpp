#include <iostream>
#include "TParticle.h"
#include <time.h>

double rnd() {
	return rand() % 1000;
}
void main() {
	srand(time(0));

	std::vector<double> tmp = {2, 2, 2};
	MyVector E;
	E = tmp;

	tmp = { 0.1, 0.1, 0.1 };
	MyVector B;
	B = tmp;

	TParticle* parts[100];
	for (int i = 0; i < 100; i++) {
		parts[i] = new TParticle(rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
	}
	//TParticle part(0, 0, 0, 1, 1, 1, 1);

	time_t start = clock();
	for (int i = 0; i < 100; i++) {
		parts[i]->updateAllAndPrint(E, B);
	}
	time_t end = clock();

	double time = ((double)end - start) / ((double)CLOCKS_PER_SEC);
	std::cout << "Execution time: " << time << " seconds.\n";
}