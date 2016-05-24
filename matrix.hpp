#ifndef HH_MATRIX_HH
#define HH_MATRIX_HH
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
#include <Eigen/Dense>

using namespace Geometry;
using namespace std;
using namespace Eigen;

class MyMatrix {

public:
	MyMatrix()=default;
	MyMatrix(int rows,int cols): rows(rows),cols(cols),M(MatrixXd(rows,cols)) {};

	MyMatrix (MyMatrix const & m)=default;

	MyMatrix & operator = (const MyMatrix & m)=default;

	double & operator() (int i, int j) {return M(i,j);};
	double operator() (int i, int j) const {return M(i,j);};

	unsigned int GetRows() const {return rows;};
	unsigned int GetCols() const {return cols;};

	friend ostream & operator<< (ostream & out, MyMatrix const & m){
		for(unsigned int i=0; i<m.GetRows(); i++){
			for (unsigned int j=0; j<m.GetCols(); j++)
				out<<m(i,j);
			out<<endl;
		}
		return out;
	};


private:
	int rows;
	int cols;
	MatrixXd M;
};

#endif