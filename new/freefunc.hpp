#ifndef _FREEFUNC_HPP_
#define _FREEFUNC_HPP_
#include <vector>
#include <string>
#include <memory>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <sstream>
#include "Mesh.hpp"
#include <fstream>
#include <algorithm>
#include <set>
#include "Geometry.hpp"

std::vector<std::array<int,2> > Polynomials(int k);
void computeDOF(std::vector<Point> const & points,unsigned int k,std::vector<double> & weights,std::vector<Point> & nodes);

#endif