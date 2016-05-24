#include "Polygon.hpp"
#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <limits>
namespace Geometry{
  
  double distance(Point2D const & a, Point2D const & b){
    // Not very efficient. This function should be implemented either
    // as friend or as a method, to access private members.
    Point2D tmp(a-b);
    return std::sqrt(tmp.x()*tmp.x()+
		     tmp.y()*tmp.y()
		     );
  }
  // ********************* BASE CLASS **********************
  
  AbstractPolygon::AbstractPolygon(std::vector<unsigned int> const & v, Vertices * p, bool check)
  :position(v),pointer(p)
  {
    if (check) this->checkConvexity();
  }
/*
  void AbstractPolygon::printCoord() const {
  	if (pointer==nullptr) {std::cout<<"Nullo"; return;}
	std::cout<<"Coordinate dei vertici:"<<std::endl;
	for(auto i : *pointer) std::cout<<i.x()<<" "<<i.y()<<std::endl; 
	}
*/
  
  void AbstractPolygon::showMe(std::ostream & out)const
  {
    if (this->size()==0){
      out<< "empty polygon" << std::endl;
    }
    else {
      out<<"Index of the vertexes:"<<std::endl;
      //out<<"    X     " << "   Y    "<<std::endl;
      for (auto const & i : this->position)
	out<<i<<" "<<std::endl;
    }
    if(this->isconvex) std::cout<<" Polygon is convex"<<std::endl;
    else std::cout<<" Polygon is not convex"<<std::endl;
  }
  
  void AbstractPolygon::checkConvexity()
  {
    std::vector<unsigned int> const & myV(this->position);
    auto mysize=this->size();
    // We consider segments and triangles as convex
    if (mysize <= 3) {
      this->isconvex=true;
      return;
    }
    //! Since we are dealing with floating points it is better to have
    //  a small number so that |a| < smallNumber means for us a==0
    double const smallNumber(1000*std::numeric_limits<double>::min());
    Point2D p;
    Point2D v;
    Point2D u;
    double res(0.0);
    double newres(0.0);
    //! C++11 sintax. decltype(expr) returns the type of the expression
    for ( decltype(mysize) i=0; i < mysize; ++i)
      {
	//p = myV[i];
      p=(*pointer)[myV[i]];
	// ! next point
	//Point2D tmp = myV[(i+1) % myV.size()];
      Point2D tmp = (*pointer)[myV[(i+1)%myV.size()]];
	v = tmp - p;
	//! next next point
	u = (*pointer)[myV[(i+2)%myV.size()]];
	if (i == 0) // in first loop direction is unknown, so save it in res
	  res = u.x() * v.y() - u.y() * v.x() + v.x() * p.y() - v.y() * p.x();
	else{
	  newres = u.x() * v.y() - u.y() * v.x() + v.x() * p.y() - v.y() * p.x();
	  if (std::abs(res)<smallNumber){
	    // The two edges are aligned, skip test and update res
	    res=newres;
	  }
	  else if ((newres > 0 && res < 0) || (newres < 0 && res > 0) )
	    {
	      this->isconvex=false;
	      return;
	    } 
	}
      }// end for
    this->isconvex=true;
    return;
  }
  
  
  using PosVector=std::vector<unsigned int>;
  // ****   POLYGON
  
  Polygon::Polygon(PosVector const & v, Vertices *p): AbstractPolygon(v,p) {}
  
  
  //! To compute the area of a polygon we use the divergence theorem.
  /*!
    \f$ int_P d\Omega = 1/2 \int{\partial P} xn_x d\gamma\f$
    The integral is computed by using trapezoidal rule on the polygon sides.
  */
  double Polygon::area() const{
    auto siz=this->size();
    if (siz<3) return 0;
    double result(0);
    PosVector const & Ver(this->position);
    // I am using C++11 sintax. decltype(expr) returns the
    // type of the expression expr. (I could have used Vertices::size_type, 
    // or just int, here I wanted to show some features of the new standard).
    for (decltype(siz) i=0; i<siz;++i){
      // Current point I use references to avoid unnecessary construction
      // of objects
      Point2D const & p1((*pointer)[Ver[i]]);
      // Other point
      //Point2D const & p2(Ver[(i+1) % siz]);
      //Point2D const & p0(Ver[(i-1) % siz]);
      Point2D const & p2((*pointer)[Ver[(i+1)%siz]]);
	  Point2D const & p0((*pointer)[Ver[(i-1)%siz]]);
      result+=p1.x()*(p2.y()-p0.y());
    }
    return 0.5*result;
  }
  
  void Polygon::showMe(std::ostream & out)const
  {
    std::cout<<" A Generic Polygon"<<std::endl;
    AbstractPolygon::showMe(out);
  }

  // ********************* SQUARE **********************
  Square::Square(PosVector const & v, Vertices *p): AbstractPolygon(v,p,false) 
  {
  	
    this->isconvex=true;
    if(v.size() != 4){
      throw std::runtime_error(" A square must be created giving four vertices");
    }

    // Check if it is a square!
    double l1=distance((*pointer)[v[1]],(*pointer)[v[0]]);
    double l2=distance((*pointer)[v[2]],(*pointer)[v[3]]);
    auto ratio=std::abs(this->area())/(l1*l2);
    if (std::abs(ratio - 1.0)>10*std::numeric_limits<double>::epsilon())
      throw std::runtime_error("Vertexes do not define a square");
     
  }
/*
  Square::Square(unsigned int origin, double length, double angle){
    this->isconvex=true;
    this->vertexes.reserve(4);
    position.push_back(origin);
    double c(std::cos(angle));
    double s(std::sin(angle));
    vertexes.push_back(Point2D(origin.x()+length*c,origin.y()+length*s));
    vertexes.push_back(Point2D(origin.x()+length*(c-s),origin.y()+length*(c+s)));
    vertexes.push_back(Point2D(origin.x()-length*s,origin.y()+length*c));
  }
*/
  
  double Square::area() const{
    if(this->size()==0) return 0.0;
    // I want the area with sign, positive if the quad is 
    // oriented counterclockwise. So the easiest thing is to use the 
    // formula using cross product (even if it is not the most efficient choice)
    Point2D v((*pointer)[position[1]]-(*pointer)[position[0]]);
    Point2D w((*pointer)[position[2]]-(*pointer)[position[0]]);
    //Point2D w(this->vertexes[2]-this->vertexes[0]);
    // area = v \times w. Positive if square counterclockwise oriented
    return v.x()*w.y()-v.y()*w.x();
    ;}
  
  void Square::showMe(std::ostream & out) const
  {
    out<<"A Square"<<std::endl;
    AbstractPolygon::showMe(out);
  }
  
  
  //********************* TRIANGLE **********************
  
  Triangle::Triangle(PosVector const & v, Vertices * p):AbstractPolygon(v,p,false){
    this->isconvex=true;    
    // Check if we give 3 vertices
    // We may use assert, in this case we would disable the control
    // in the released version (-DNDEBUG). We prefer here to exit the program
    if(this->size() != 3){
      throw std::runtime_error(" A triangle must be created giving three vertices");
    }
  }
  

  double Triangle::area() const
  {
    if(this->size()==0) return 0.0;
    // I use the cross product since this is a 
    // signed area!
    Point2D v((*pointer)[position[1]]-(*pointer)[position[0]]);
    Point2D w((*pointer)[position[2]]-(*pointer)[position[0]]);
    // area = 0.5* v \times w. Positive if triangle counterclockwise oriented
    return 0.5*(v.x()*w.y()-v.y()*w.x());
    ;}
  
  void Triangle::showMe(std::ostream & out) const{
    out<<"A Triangle"<<std::endl;
    AbstractPolygon::showMe(out);
  }
}
