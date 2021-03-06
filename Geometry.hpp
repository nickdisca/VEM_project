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

	//get the coordinates
	void getCoordinates(double & x, double & y) const {
		x=M_coor[0]; y=M_coor[1]; }

	//operator []
	double const operator[] (int i) const {return M_coor[i];}
	double & operator[] (int i) {return M_coor[i];}

	//operations
	friend Point operator +(const Point &, const Point &);
	friend Point operator -(const Point &, const Point &);
	Point operator *(const double &)const;
	friend Point operator*(const double &, const Point &);

	//comparison operator between two points
	friend bool operator == (Point const & a, Point const & b){return ((a[0]==b[0]) && (a[1]==b[1]));};
	//friend bool operator == (Point const & a, Point const & b){return (std::abs(a[0]-b[0])<=1e-8 && std::abs(a[1]-b[1])<=1e-8);};
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
	
	//compute normal with respect to edge i
	Point Normal(unsigned int edge_num);

	//given two vectors of uint and points, assigns them to the class members
	void setDof(std::vector<unsigned int> const &, std::vector<Point> *);
	//get the coordinates of the dofs (GLL points)
	std::vector<Point> getDof() const;
	//get the indexes of the dofs (GLL points)
	std::vector<unsigned int> getBDindexes(){return dof;};

	//compute local matrices
	MatrixType ComputeD(unsigned int k);
	MatrixType ComputeB(unsigned int k);
	MatrixType ComputeG(unsigned int k);
	MatrixType LocalStiffness(unsigned int k);
	MatrixType ComputeH(unsigned int k, std::function<double (double,double)> weight, unsigned int krows, unsigned int kcols);
	MatrixType ComputeH(unsigned int k, std::function<double (double,double)> weight=[](double x, double y) {return 1.0;}) 
		{return ComputeH(k,weight,k,k);}	
	MatrixType ComputeC(unsigned int k);
	MatrixType LoadTerm(unsigned int k,std::function<double (double,double)> f);
	MatrixType LocalMass(unsigned int k);
	MatrixType LocalStiffness_weighted
		(unsigned int k, std::function<double (double,double)> mu, double mu_bar, bool constant_mu);
	MatrixType LocalTransport
		(unsigned int k, std::function<double (double,double)> beta_x,std::function<double (double,double)> beta_y);
	MatrixType ComputeE(unsigned int k, unsigned int VAR);

	//compute the dofs for a function uex (i.e. compute its VEM approximation)
	MatrixType LocalConvert(unsigned int k, std::function<double (double,double)> uex);

	//output
	friend std::ostream & operator << (std::ostream &, const Polygon &);

private:
	std::vector<unsigned int> vertexes;
	std::vector<Point> * pointer;
	std::vector<unsigned int> dof;
	std::vector<Point> * pointer_dof;
};

#endif
