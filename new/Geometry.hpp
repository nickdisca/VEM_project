#ifndef __GEOMETRY_HPP_
#define __GEOMETRY_HPP_
#include <array>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

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

	//confronto
	friend bool operator == (Point const & a, Point const & b){return ((a[0]==b[0]) && (a[1]==b[1]));};
	friend bool operator <(Point const &f, Point const &s){
		if (f[0]==s[0]) return f[1]<s[1];
		return f[0]<s[0];
	};

	//distance between two points
	friend double distance(const Point &, const Point &);

	//output
	friend std::ostream & operator << (std::ostream &, const Point &);

private:
	std::array<double,2> M_coor;
};


inline  Point operator - (const Point & a, const Point & b)
	{ return Point(a.M_coor[0]-b.M_coor[0],a.M_coor[1]-b.M_coor[1]); };

inline Point operator + (const Point & a, const Point & b)
	{ return Point(a.M_coor[0]+b.M_coor[0],a.M_coor[1]+b.M_coor[1]); };


inline  Point operator *(const double & d, const Point & p)
	{ return p*d; };

//END OF CLASS POINT


using MatrixType=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;

class Polygon {
public:
	//standard methods
	Polygon():vertexes(),pointer(nullptr){};
	Polygon(const Polygon &)=default; 
	Polygon & operator=(const Polygon&)=default; 
	Polygon(Polygon &&)=default;
	Polygon & operator=(Polygon&&)=default;
	//~Polygon(){};

	//constructor giving the two vectors
	Polygon(std::vector<unsigned int> const v, std::vector<Point> * p):vertexes(v),pointer(p){};

	//return the indexes
	std::vector<unsigned int> getVertexes(){return vertexes;};
	//get the points
	std::vector<Point> getPoints() const;
	//size of the polygon
	unsigned int const size() const {return vertexes.size();};

	//operator []
	unsigned int const operator[] (int i) const {return vertexes[i];}
	unsigned int & operator[] (int i) {return vertexes[i];}

	//diameter, area, centroid
	double area() const;
	Point centroid() const;
	double diameter() const;
	
	Point Normal(unsigned int edge_num);

	//boundary dof
	void setDof(std::vector<unsigned int> const &, std::vector<Point> *);
	std::vector<Point> getDof() const;

	//local matrices
	MatrixType ComputeD(unsigned int k);
	MatrixType ComputeB(unsigned int k);
	MatrixType ComputeG(unsigned int k);

	//output
	friend std::ostream & operator << (std::ostream &, const Polygon &);

private:
	std::vector<unsigned int> vertexes;
	std::vector<Point> * pointer;
	std::vector<unsigned int> dof;
	std::vector<Point> * pointer_dof;
};

#endif
