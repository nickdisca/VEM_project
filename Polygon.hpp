#ifndef HH_POLYGON_HH
#define HH_POLYGON_HH
#include "matrix.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <utility>

using MatrixType=Matrix<double,Dynamic,Dynamic>;

namespace Geometry
{
  class Point2D
  {
  public:
    //! Constructor giving coordinates.
    Point2D(double xx=0.0, double yy=0.0):coor{{xx,yy}}{}
    //! Copy constructor
    Point2D(const Point2D&)=default;
    //! Returns coordinates in a array<double>.
    std::array<double,2> get() const { return coor;}
    //! Sets point coordinates
    void set(double const &xx, double const &yy)
    {
      coor={{xx,yy}};
    }
    //! x coordinate
    double x() const {return coor[0];}
    //! y coordinate
    double y() const {return coor[1];}
    //! Subtraction is implemented as an external friend function
    friend Point2D operator - (Point2D const & a, Point2D const & b);
    //! Addition is implemented as an external friend function
    friend Point2D operator + (Point2D const & a, Point2D const & b);

    friend bool operator == (Point2D const & a, Point2D const & b){return ((a.x()==b.x()) && (a.y()==b.y()));};
  private:
    std::array<double,2> coor;
  };

  //! An alias
  using R2Vector=Point2D;

  //! subtraction operator (inline)
  inline Point2D operator - (Point2D const & a, Point2D const & b){
    return Point2D(a.coor[0]-b.coor[0],
                   a.coor[1]-b.coor[1]);
  }

  //! Addition operator (inline)
  inline Point2D operator + (Point2D const & a, Point2D const & b){
    return Point2D(a.coor[0]+b.coor[0],
                   a.coor[1]+b.coor[1]);
  }

  //! Distance between points
  double distance(Point2D const & a, Point2D const & b);
 
  //! Polygon vertices are just vectors of points.
  using Vertices=std::vector<Point2D>;

  //! Defines the common interface of polygons.
  class AbstractPolygon
  {
  public:
    //! Constructor taking vertices
    /*! 
      It checks convexity if check=true
     */
    AbstractPolygon(Vertices const & v, bool check=true);
    //! Default constructor is defaulted
    /*! 
      It is up to the derived classes to fill the vertexex and other info correctly
    */
    AbstractPolygon()=default;
    //! Assignment
    AbstractPolygon & operator=(AbstractPolygon const&)=default;
    //! Copy constructor
    AbstractPolygon(AbstractPolygon const &)=default;
    //! Move constructor
    AbstractPolygon(AbstractPolygon&&)=default;
    //! Move constructor
    AbstractPolygon & operator=(AbstractPolygon&&)=default;
    //! virtual destructor
    virtual ~AbstractPolygon(){};
    //! Returns the number of vertices.
    /*!  We return Vertices::size_type and not just int because
      size_type is guaranteed to be the correct type for indexes in
      stl vectors. Its actual type may be implementation dependent.
      In this case, however, int would have been fine (size_type is
      guaranteed to be an integral type, more precisely
      a type convertible to unsigned int).
    */
    Vertices::size_type size() const {return vertexes.size();}
    //! Is the polygon convex?
    bool isConvex() const {return isconvex;}
    //! Returns the vertices (read only)
    Vertices const & theVertices()const {return vertexes;}
    //! Outputs some info on the polygon
    virtual void showMe(std::ostream & out=std::cout) const;
    //! The area of the polygon (with sign!).
    /*!
      It is a pure virtual function.
      The implementation is left to the derived classes.
    */
    virtual double area() const=0;
    //! The centroid of the polygon
    Point2D Centroid() const;
    //! The diameter of the polygon
    double Diameter() const;

    Vertices BoundaryDof(int k) const;
    //std::vector<std::array<int,2> > Polynomials(int k) const;
    MatrixType ComputeD(int k);
    MatrixType ComputeB(int k);
    MatrixType ComputeG(int k);
    MatrixType ComputeGTilda(int k);
    MatrixType ComputeStiffness(int k);
    double ComputeIntegral(int k, int d1, int d2);
    double IntegralWithDof(int k, int edge, int phi);
    Point2D Normal(int edge);

  protected:
    Vertices vertexes;
    bool isconvex;
    //! Test convexity of the polygon
    void checkConvexity();
    
  };

  //! Class for a generic Polygon
  /*
    A generic Polygon is defined by a set of Vertices which are
    provided by the user
   */
  class Polygon: public AbstractPolygon
  {
  public:
    //! Default constructor.
    //! Polygon may be constructed giving Vertices;
    Polygon(Vertices const & v);
    //! Destructor
    virtual ~Polygon(){};
    /*!
      The area is positive if vertices are given in 
      counterclockwise order
    */
    virtual double area() const;
    //! Specialised version for generic polygons.
    virtual void showMe(std::ostream & out=std::cout) const;

  };

  //! A square
  /*!
    The square is a final class derived from polygon.
   */
  class Square final: public AbstractPolygon
  {
  public:
    Square(Vertices const & v);
    //!Special constructor valid only for squares.
    /*!
      /param origin Point which gives the first vertex of the square.
      /param length The length of the side.
      /param angle In radians, tells how the square is  rotated. 
     */
    Square(Point2D origin, double length,double angle=0.0);
    Square(Square const &)=default;
    Square(Square&&)=default;
    Square & operator=(const Square &)=default;
    Square & operator=(Square &&)=default;
    //! Specialised version for squares
    double area() const;
    //! Specialised version for squares.
    void showMe(std::ostream & out=std::cout) const;
  };
  
  //! A triangle
  class Triangle final: public AbstractPolygon
  {
  public:
    Triangle(Vertices const &);
    Triangle(Triangle const &)=default;
    Triangle(Triangle&&)=default;
    Triangle & operator=(const Triangle &)=default;
    Triangle & operator=(Triangle &&)=default;
    //! Specialised for Triangles
    double area() const;
    //! Specialised for Triangles
    virtual void showMe(std::ostream & out=std::cout) const;
  };
  
}

#endif
