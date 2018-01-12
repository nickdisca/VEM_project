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
	Quadrature()=default;
	Quadrature(const Polygon & PP): P(PP){};
	Quadrature(const Quadrature &)=default; 
	Quadrature & operator=(const Quadrature&)=default; 
	Quadrature(Quadrature &&)=default;
	Quadrature & operator=(Quadrature&&)=default;
	~Quadrature(){};

	//output
	friend std::ostream & operator << (std::ostream &, const Quadrature &);

	//divide in triangles (restituisce un vettore di "triangoli")
	std::vector<std::array<Point,3> > divide();

	//calcola integrale globale sul poligono
	double local_int(std::array<Point,3> & tria, std::function<double(double,double)> f, unsigned int n);
	double global_int(std::function<double(double,double)> f, unsigned int n);

	//mappa e valutazione
	double map(std::array<Point,3> & p, MatrixType & B, MatrixType & b);
	Point map_eval(Point & x,MatrixType & B, MatrixType & b);


private:
	Polygon P;
};

void Gauss_Leg(double a, double b, unsigned int n, std::vector<double> & nodes,std::vector<double> & weights);
void reference(unsigned int n,
	std::vector<double> & nodes1d, std::vector<double> & weights1d, std::vector<Point> & nodes2d, std::vector<double> & weights2d);

#endif