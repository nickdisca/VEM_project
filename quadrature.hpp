#ifndef __QUADRATURE_HPP_
#define __QUADRATURE_HPP_
#include <array>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <functional>
#include "Geometry.hpp"

class Quadrature {
public:

	//standard
	Quadrature()=default;
	Quadrature(const Polygon & PP): P(PP){};
	Quadrature(const Quadrature &)=default; 
	Quadrature & operator=(const Quadrature&)=default; 
	Quadrature(Quadrature &&)=default;
	Quadrature & operator=(Quadrature&&)=default;
	~Quadrature(){};

	//output
	friend std::ostream & operator << (std::ostream &, const Quadrature &);

	//divide in triangles (returns a vector of "triangles")
	std::vector<std::array<Point,3> > divide();

	//compute integral on a triangle
	double local_int(std::array<Point,3> & tria, std::function<double(double,double)> f, unsigned int n);
	//compute integral on the polygon
	double global_int(std::function<double(double,double)> f, unsigned int n);

	//compute affine map and evaluate it (matrix-vector multiplication)
	double map(std::array<Point,3> & p, MatrixType & B, MatrixType & b);
	Point map_eval(Point & x,MatrixType & B, MatrixType & b);

private:
	Polygon P;
};

//compute 1D GL nodes and weights on a generic interval [a,b]
void Gauss_Leg(double a, double b, unsigned int n, std::vector<double> & nodes,std::vector<double> & weights);
//compute 1D, 2D GL nodes and weights on reference element
void reference(unsigned int n,
	std::vector<double> & nodes1d, std::vector<double> & weights1d, std::vector<Point> & nodes2d, std::vector<double> & weights2d);

#endif
