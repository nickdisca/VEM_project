#ifndef __GEOMETRY_HPP_
#define __GEOMETRY_HPP_
#include <array>
#include <iostream>

//class representing a point in 2D
class Point {
public:
	//standard methods
	explicit  Point(double x=0.0, double y=0.0):M_coor{{x,y}}{};
	Point(const Point &)=default; 
	Point & operator=(const Point&)=default; 
	Point(Point &&)=default;
	Point & operator=(Point&&)=default;
	~Point(){};
	
	//set the coordinates
	void setCoordinates(double x, double y){
		M_coor[0]=x; M_coor[1]=y;}

	//get the coordinates (stored in the input elements, not as return value)
	void getCoordinates(double & x, double & y) const {
		x=M_coor[0]; y=M_coor[1]; }

	//operator [], const and not const version
	double const operator[] (int i) const {return M_coor[i];}
	double & operator[] (int i) {return M_coor[i];}

	//operations
	friend Point operator +(const Point &, const Point &);
	friend Point operator -(const Point &, const Point &);
	Point operator *(const double &)const;
	friend Point operator*(const double &, const Point &);

	//distance between two points
	friend double distance(const Point &, const Point &);

	//output
	friend std::ostream & operator << (std::ostream &, const Point &);

protected:
	std::array<double,2> M_coor;
};


inline  Point operator - (const Point & a, const Point & b)
	{ return Point(a.M_coor[0]-b.M_coor[0],a.M_coor[1]-b.M_coor[1]); };

inline Point operator + (const Point & a, const Point & b)
	{ return Point(a.M_coor[0]+b.M_coor[0],a.M_coor[1]+b.M_coor[1]); };


inline  Point operator *(const double & d, const Point & p)
	{ return p*d; };

//END OF CLASS POINT



#endif
