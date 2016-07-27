#ifndef __QUADRATURE_HPP_
#define __QUADRATURE_HPP_
#include <array>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "Geometry.hpp"

class Quadrature {
public:
	Quadrature()=default;
	Quadrature(const Polygon & PP): P(PP){};
	Quadrature(const Quadrature &)=default; 
	Quadrature & operator=(const Quadrature&)=default; 
	Quadrature(Quadrature &&)=default;
	Quadrature & operator=(Quadrature&&)=default;
	~Quadrature(){};

	//output
	friend std::ostream & operator << (std::ostream &, const Quadrature &);


private:
	Polygon P;
};

void Gauss_Leg(double a, double b, unsigned int n, std::vector<double> & nodes,std::vector<double> & weights);
void reference(unsigned int n,
	std::vector<double> & nodes1d, std::vector<double> & weights1d, std::vector<Point> & nodes2d, std::vector<double> & weights2d);

#endif