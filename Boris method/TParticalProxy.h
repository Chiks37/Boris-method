#include <iostream>
#include <vector>
#pragma once
class TParticalProxy
{
	//Необходимость использование double или float требует дальнейшего исследования
	//позиция и момент - вектора трёх элементов!
	std::vector<double>* position; // Particle position (x, y, z)
	std::vector<double>* momentum; // Particle momentum (px, py, pz)
	double* weight; // Particle weight
	double* gamma; // Particle γ-factor
	short int* type; // Particle type
};

