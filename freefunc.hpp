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

//computes the discretization degree in the correct order x^{s_1}*y^{s_2}
std::vector<std::array<int,2> > Polynomials(int k);

//computes the degrees of freedom (i.e. GLL points) on edges)
void computeDOF(std::vector<Point> const & points,unsigned int k,std::vector<double> & weights,std::vector<Point> & nodes);

#endif
