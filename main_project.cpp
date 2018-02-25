#include <iostream>
#include <vector>
#include <fstream>
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "quadrature.hpp"
#include "freefunc.hpp"

using namespace std;

int main()
{

//set polynomial degree
int k=1;

//read mesh from file structured as [list of 2D points; connectivity matrix; indexes of boundary points]
std::string str="./Meshes/Uniform/64.dat";
MeshReader read(false);
Mesh m(str,read,k);
//cout<<m;
cout<<"Total area = "<<m.area()<<endl;
cout<<"Maximum diameter = "<<m.max_diam()<<endl;

//example
//exact solution: sin(2*pi*x)*sin(2*pi*y) for elliptic problem
double PI=4.0*std::atan(1.0);
auto f=[PI](double x,double y){return -2.0*PI*(std::cos(2.0*PI*x)*std::sin(2.0*PI*y)-
	2.0*PI*(x+1.0)*std::sin(2.0*PI*x)*std::sin(2.0*PI*y)-2*PI*(x+1.0)*std::sin(2.0*PI*x)*std::sin(2.0*PI*y))+2.0*PI*
	(std::cos(2.0*PI*x)*std::sin(2.0*PI*y)+std::sin(2.0*PI*x)*std::cos(2.0*PI*y));};
auto g=[PI](double x,double y){return std::sin(2.0*PI*x)*std::sin(2.0*PI*y);};
auto uex=[PI](double x,double y){return std::sin(2.0*PI*x)*std::sin(2.0*PI*y);};
auto mu=[](double x, double y) {return x+1.0;}; double mu_bar=1.5;
auto beta_x=[](double x, double y) {return 1.0;};
auto beta_y=[](double x, double y) {return 1.0;};


//solve the problem
MatrixType U=m.solve(f,g, mu, mu_bar, false, beta_x, beta_y);

//compute the norms
MatrixType UEX=m.VEMConvert(uex);
m.Allnorms(U,UEX);

//output
//VEM solution is stored in the file output.dat as [pointX pointY value(point)]

return 0;
}
