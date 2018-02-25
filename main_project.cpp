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
std::string str="./mesh_example.dat";
MeshReader read(false);
Mesh m(str,read,k);
//cout<<m;
cout<<"Total area = "<<m.area()<<endl;
cout<<"Maximum diameter = "<<m.max_diam()<<endl;

//exact solution: x+y for Laplace problem
auto f=[](double x,double y){return 0.0;};
auto g=[](double x,double y){return x+y;};
auto uex=[](double x,double y){return x+y;};
auto mu=[](double x, double y) {return 1.0;}; double mu_bar=1.0;
auto beta_x=[](double x, double y) {return 0.;};
auto beta_y=[](double x, double y) {return 0.;};


//solve the problem: the bool true if the diffusion coefficient is constant, false otherwise
MatrixType U=m.solve(f,g, mu, mu_bar, true, beta_x, beta_y);

//compute the norms
MatrixType UEX=m.VEMConvert(uex);
m.Allnorms(U,UEX);

//output
//VEM solution is stored in the file output.dat as [pointX pointY value(point)]

return 0;
}
