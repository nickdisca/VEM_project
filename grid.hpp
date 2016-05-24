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
#include <algorithm>

using namespace Geometry;
using namespace std;

class Grid
  {
  public:
    Grid ()=default;
    Grid (const Grid & )=default;
    Grid & operator= (const Grid & )= default;
    Grid (ifstream & file);

    double area();
    unsigned int grid_size(){return abspol.size();};
    void printedges();
    void printedgesIndex(ofstream & ost1, ofstream & ost2, ofstream & ost3);

  private:
    std::vector<Point2D> vect;
    std::vector<std::shared_ptr<AbstractPolygon> > abspol;
    std::vector<Edge> AllEdges,Boundary;
  };

  #endif