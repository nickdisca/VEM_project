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


std::string str="Meshes/Non_uniform/128.dat";
MeshReader read(false);
Mesh m(str,read,2);
cout<<m;
cout<<"Total area = "<<m.area()<<endl;
cout<<"Maximum diameter = "<<m.max_diam()<<endl;
//m.boundaryDOF();
//cout<<m;
//MatrixType K=m.GlobalStiffness();
//K=m.GlobalStiffness();
//cout<<K;

/*
//test with unit square (k=2)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(1.0,0.0)); poi.push_back(Point(1.0,1.0)); poi.push_back(Point(0.0,1.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); 
Polygon p(line,&poi);
std::vector<Point> dof; 
dof.push_back(Point(0.5,0.0)); dof.push_back(Point(1.0,0.5)); dof.push_back(Point(0.5,1.0)); dof.push_back(Point(0.0,0.5));
p.setDof(line,&dof);
cout<<p;
cout<<"Matrix B:"<<endl<<p.ComputeB(2)<<endl;
//cout<<"Matrix D:"<<endl<<p.ComputeD(2)<<endl;
//cout<<"Matrix G:"<<endl<<p.ComputeG(2)<<endl;
//cout<<"Matrix K (local stiffness):"<<endl<<p. LocalStiffness(2)<<endl;
//cout<<"Matrix H:"<<endl<<p.ComputeH(2)<<endl;
//cout<<"Matrix C:"<<endl<<p.ComputeC(2)<<endl;
//auto ff=[](double x,double y){return 1.0;};
//cout<<"Vector F:"<<endl<<p.LoadTerm(2,ff)<<endl;
//cout<<(p.ComputeH(2).lu()).solve(p.ComputeC(2));
*/

/*
//test with unit square (k=3)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(1.0,0.0)); poi.push_back(Point(1.0,1.0)); poi.push_back(Point(0.0,1.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); 
Polygon p(line,&poi);
std::vector<Point> dof; 
double radice=1./std::sqrt(5);
dof.push_back(Point(0.5*(1-radice),0.0)); dof.push_back(Point(0.5*(1+radice),0.0));
dof.push_back(Point(1.0,0.5*(1-radice))); dof.push_back(Point(1.0,0.5*(1+radice)));
dof.push_back(Point(0.5*(1+radice),1.0)); dof.push_back(Point(0.5*(1-radice),1.0));
dof.push_back(Point(0.0,0.5*(1+radice))); dof.push_back(Point(0.0,0.5*(1-radice)));
line.push_back(4); line.push_back(5); line.push_back(6); line.push_back(7); 
p.setDof(line,&dof);
cout<<p;
MatrixType B=p.ComputeB(3), D=p.ComputeD(3), G=p.ComputeG(3);
//cout<<"Matrix B:"<<endl<<p.ComputeB(3)<<endl;
//cout<<"Matrix D:"<<endl<<p.ComputeD(3)<<endl;
//MatrixType B=p.ComputeB(3); MatrixType D=p.ComputeD(3); MatrixType M=B*D;
//cout<<"Matrix product:"<<endl<<B*D<<endl;
//cout<<"Matrix product:"<<endl<<MM<<endl;
cout<<"Matrix K (local stiffness):"<<endl<<p.LocalStiffness(3)<<endl;
cout<<"Matrix K (local stiffness):"<<endl<<p.LocalStiffness(3)<<endl;
///B=p.ComputeB(3); D=p.ComputeD(3); G=p.ComputeG(3);
//cout<<"Matrix B:"<<endl<<p.ComputeB(3)<<endl;
//cout<<"Matrix D:"<<endl<<p.ComputeD(3)<<endl;
//B=p.ComputeB(3); D=p.ComputeD(3); MatrixType GG=p.ComputeB(3)*p.ComputeD(3);
cout<<"Matrix B:"<<endl<<p.ComputeB(3)<<endl;
cout<<"Matrix D:"<<endl<<p.ComputeD(3)<<endl;;
B=p.ComputeB(3); D=p.ComputeD(3);
cout<<"Matrix B:"<<endl<<B<<endl;
cout<<"Matrix D:"<<endl<<D<<endl;
cout<<"Matrix G:"<<endl<<B*D<<endl;

//cout<<"Matrix H:"<<endl<<p.ComputeH(3)<<endl;
//cout<<"Matrix C:"<<endl<<p.ComputeC(3)<<endl;
//auto f=[](double x,double y){return 1.0;};
//cout<<"Vector F:"<<endl<<p.LoadTerm(3,f)<<endl;
//cout<<(p.ComputeH(2).lu()).solve(p.ComputeC(2));
*/

//auto f=[](double x,double y){return 0.0;};
//auto g=[](double x,double y){return 1.0;};

//double pi=4.0*std::atan(1.0);
//auto f=[pi](double x,double y){return 15.0*std::sin(pi*x)*std::sin(pi*y);};
//auto g=[pi](double x,double y){return (1.0-x)*y*std::sin(pi*x);};

//exact solution: 1
//auto f=[](double x,double y){return 0.0;};
//auto g=[](double x,double y){return 1.0;};
//auto uex=[](double x,double y){return 1.0;};

//exact solution: x^2+y^2
//auto f=[](double x,double y){return -4.0;};
//auto g=[](double x,double y){return x*x+y*y;};
//auto uex=[](double x,double y){return x*x+y*y;};

//exact solution: sin(x)*sin(y)
//auto f=[](double x,double y){return 2.0*std::sin(x)*std::sin(y);};
//auto g=[](double x,double y){return std::sin(x)*std::sin(y);};
//auto uex=[](double x,double y){return std::sin(x)*std::sin(y);};

//exact solution: x^2*abs(x) centered in x=1/2
auto f=[](double x,double y){return -6.0*std::abs(x-1./2);};
auto g=[](double x,double y){return (x-1./2)*(x-1./2)*std::abs(x-1./2);};
auto uex=[](double x,double y){return (x-1./2)*(x-1./2)*std::abs(x-1./2);};

//exact solution: x*abs(x) centered in x=1/2
//auto f=[](double x,double y){return -2.0*((x>1./2)?1.0:(x < 1./2)?-1.0:0.);};
//auto g=[](double x,double y){return (x-1./2)*std::abs(x-1./2);};
//auto uex=[](double x,double y){return (x-1./2)*std::abs(x-1./2);};


//MatrixType F=m.GlobalLoad(f);
//cout<<F<<endl;
MatrixType U=m.solve(f,g);
//cout<<"Solution"<<U<<endl;

MatrixType UEX=m.VEMConvert(uex);
//cout<<"Vem approximation of exact solution"<<endl;
//cout<<UEX<<endl;
m.Allnorms(U,UEX);
//cout<<m;
//m.Dirichlet();



/*
//test with unit square (k=2)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(1.0,0.0)); poi.push_back(Point(1.0,1.0)); poi.push_back(Point(0.0,1.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); 
Polygon p(line,&poi);
std::vector<Point> dof; 
dof.push_back(Point(0.5,0.0)); dof.push_back(Point(1.0,0.5)); dof.push_back(Point(0.5,1.0)); dof.push_back(Point(0.0,0.5));
p.setDof(line,&dof);
cout<<p;
cout<<"Matrix B:"<<endl<<p.ComputeB(2)<<endl;
cout<<"Matrix D:"<<endl<<p.ComputeD(2)<<endl;
cout<<"Matrix G:"<<endl<<p.ComputeG(2)<<endl;
cout<<"Matrix K (local stiffness):"<<endl<<p. LocalStiffness(2)<<endl;
cout<<"Matrix H:"<<endl<<p.ComputeH(2)<<endl;
cout<<"Matrix C:"<<endl<<p.ComputeC(2)<<endl;
auto f=[](double x,double y){return 1.0;};
cout<<"Vector F:"<<endl<<p.LoadTerm(2,f)<<endl;
//cout<<(p.ComputeH(2).lu()).solve(p.ComputeC(2));
*/
/*
//test with unit square (k=1)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(1.0,0.0)); poi.push_back(Point(1.0,1.0)); poi.push_back(Point(0.0,1.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); 
Polygon p(line,&poi);
std::vector<Point> dof; 
std::vector<unsigned int> ll;
p.setDof(ll,&dof);
cout<<p;
cout<<"Matrix B:"<<endl<<p.ComputeB(1)<<endl;
cout<<"Matrix D:"<<endl<<p.ComputeD(1)<<endl;
cout<<"Matrix G:"<<endl<<p.ComputeG(1)<<endl;
cout<<"Matrix K (local stiffness):"<<endl<<p. LocalStiffness(1)<<endl;
*/

/*
//test with pentagon (k=1)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(3.0,0.0)); poi.push_back(Point(3.0,2.0)); poi.push_back(Point(1.5,4.0));
	poi.push_back(Point(0.0,4.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); line.push_back(4);
Polygon p(line,&poi);
std::vector<Point> dof; 
std::vector<unsigned int> ll;
p.setDof(ll,&dof);
cout<<p;
cout<<"Matrix B:"<<endl<<p.ComputeB(1)<<endl;
cout<<"Matrix D:"<<endl<<p.ComputeD(1)<<endl;
cout<<"Matrix G:"<<endl<<p.ComputeG(1)<<endl;
*/
/*
//test with unit pentagon (k=2)
std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(3.0,0.0)); poi.push_back(Point(3.0,2.0)); poi.push_back(Point(1.5,4.0));
	poi.push_back(Point(0.0,4.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); line.push_back(4);
Polygon p(line,&poi);
std::vector<Point> dof; 
dof.push_back(Point(1.5,0.0)); dof.push_back(Point(3.0,1.0)); dof.push_back(Point(2.25,3.0)); dof.push_back(Point(0.75,4.0));
	dof.push_back(Point(0.0,2.0));
p.setDof(line,&dof);
cout<<p;
cout<<"Matrix B:"<<endl<<p.ComputeB(2)<<endl;
cout<<"Matrix D:"<<endl<<p.ComputeD(2)<<endl;
cout<<"Matrix G:"<<endl<<p.ComputeG(2)<<endl;
*/

/*
Quadrature Q(p);
cout<<Q;

std::vector<double> nodes, weights;
Gauss_Leg(1.0,2.0,4,nodes,weights);
for (auto i : nodes) cout<<"nodi "<<i<<endl;
for (auto i : weights) cout<<"pesi "<<i<<endl;

std::vector<double> nodes1d, weights1d,weights2d;
std::vector<Point> nodes2d;
reference(3,nodes1d,weights1d,nodes2d,weights2d);
//for (auto i : nodes1d) cout<<"nodi "<<i<<endl;
//for (auto i : weights1d) cout<<"pesi "<<i<<endl;
//for (auto i : nodes2d) cout<<"nodi "<<i<<endl;
//for (auto i : weights2d) cout<<"pesi "<<i<<endl;
cout<<"The polygon is divided in these triangles:"<<endl;
auto V=Q.divide();
//auto f=[](double x,double y){return 0.5*x*y;};
auto f=[](double x,double y){return x*x;};
unsigned int n=3; //number of points on each edge
cout<<"The global integral is = "<<Q.global_int(f,n)<<endl;
*/
return 0;
}
