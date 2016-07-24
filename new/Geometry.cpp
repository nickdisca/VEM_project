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
	ost<<"X = "<<p.M_coor[0]<<" Y = "<<p.M_coor[1]<<std::endl;
	return ost;
};




std::vector<Point> Polygon::getPoints() const{
	std::vector<Point> tmp(this->size());
	for (unsigned int i=0; i<tmp.size(); ++i) 
		tmp[i]=pointer->operator[](vertexes[i]);
	return tmp;
};


std::ostream & operator << (std::ostream & ost, const Polygon & p){
	ost<<"The polygon is the following: "<<std::endl;
	for (unsigned int i=0; i<p.size(); ++i) 
		ost<<"Index "<<" = "<<p.vertexes[i]<<" corresponding to point "<<p.pointer->operator[](p.vertexes[i]);
return ost;
};


double Polygon::area() const{
	auto size=this->size();
	if (size<3) return 0.0;
	double result(0);
	auto ver(this->getPoints());

	for (decltype(size) i=0; i<size;++i){
		Point const & p1(ver[i]);
		Point const & p2(ver[(i+1) % size]);
		Point const & p0(ver[(i-1) % size]);
		result+=p1[0]*(p2[1]-p0[1]);
	}
	return 0.5*result;
};

Point Polygon::centroid() const{
	double x=0.0,y=0.0;
	auto ver(this->getPoints());
	unsigned int size=ver.size();
	for (unsigned int i=0; i<size; ++i){
	  x+=(ver[i][0]+ver[(i+1)%size][0])*(ver[i][0]*ver[(i+1)%size][1]-ver[(i+1)%size][0]*ver[i][1]);
	  y+=(ver[i][1]+ver[(i+1)%size][1])*(ver[i][0]*ver[(i+1)%size][1]-ver[(i+1)%size][0]*ver[i][1]);
	}
	x=x/(6*area());
	y=y/(6*area());
	return Point{x,y};
};

double Polygon::diameter() const{
	double d{0.0};
	auto ver(this->getPoints());
	for (unsigned int i=0; i<this->size()-1; i++){
		for (unsigned int j=0; j<this->size(); j++)
			d=std::max(d,distance(ver[i],ver[j]));
	}
	return d;
};