#ifndef HH_GRID_HH
#define HH_GRID_HH
#include "Polygon.hpp"
#include "edge.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <fstream>
#include <sstream>
#include <set>
#include <Eigen/Dense>
#include <algorithm>

using namespace Geometry;
using namespace std;
using namespace Eigen;

class Grid
  {
  public:
    Grid ()=default;
    Grid (const Grid & )=default;
    Grid & operator= (const Grid & )= default;
    Grid (ifstream & file);

    double area();
    unsigned int grid_size(){return abspol.size();};
    unsigned int boundary_size() {return boundary.size();};
    unsigned int vertices_size() {return coord.size();};

    using ReturnMat = Matrix<double,Dynamic,Dynamic>;
    //ReturnMat total_D
    

  private:
    std::vector<Point2D> coord;
    std::vector<std::shared_ptr<AbstractPolygon> > abspol;
    std::vector<unsigned int> boundary;
  };

  #endif