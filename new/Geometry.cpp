#include <cmath>
#include <iostream>
#include "Geometry.hpp"

Point Point::operator *(const double & d) const
  {return Point(d*M_coor[0],d*M_coor[1]);};

double distance(const Point & a, const Point & b){
	Point tmp(b-a);
	return std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
};

std::ostream & operator << (std::ostream & ost, const Point & p){
	ost<<"The point has coordinates X= "<<p.M_coor[0]<<" Y= "<<p.M_coor[1]<<std::endl;
	return ost;
};
