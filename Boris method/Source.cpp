#include <iostream>
#include "TParticle.h"

void main() {
	std::vector<double> E = {2, 2, 2};
	std::vector<double> B = { 0.1, 0.1, 0.1 };
	TParticle part(0, 0, 0, 1, 1, 1, 1);
	part.updateAllAndPrint(E, B);
}