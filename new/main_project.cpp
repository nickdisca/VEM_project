#include <iostream>
#include <vector>
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "quadrature.hpp"

using namespace std;

int main()
{

/*
std::string str="single.dat";
MeshReader read(false);
Mesh m(str,read,2);
//cout<<m;
cout<<"Total area = "<<m.area()<<endl;
m.boundaryDOF();
cout<<m;
*/

std::vector<Point> poi; 
poi.push_back(Point(0.0,0.0)); poi.push_back(Point(1.0,0.0)); poi.push_back(Point(1.0,1.0)); poi.push_back(Point(0.0,1.0));
std::vector<unsigned int> line;
line.push_back(0); line.push_back(1); line.push_back(2); line.push_back(3); 
Polygon p(line,&poi);
std::vector<Point> dof; 
dof.push_back(Point(0.5,0.0)); dof.push_back(Point(1.0,0.5)); dof.push_back(Point(0.5,1.0)); dof.push_back(Point(0.0,0.5));
p.setDof(line,&dof);
/*
cout<<p;
cout<<p.ComputeB(2)<<endl;
cout<<p.ComputeD(2)<<endl;
*/

Quadrature Q(p);
cout<<Q;
/*
std::vector<double> nodes, weights;
Gauss_Leg(1.0,2.0,4,nodes,weights);
for (auto i : nodes) cout<<"nodi "<<i<<endl;
for (auto i : weights) cout<<"pesi "<<i<<endl;
*/
std::vector<double> nodes1d, weights1d,weights2d;
std::vector<Point> nodes2d;
reference(3,nodes1d,weights1d,nodes2d,weights2d);
//for (auto i : nodes1d) cout<<"nodi "<<i<<endl;
//for (auto i : weights1d) cout<<"pesi "<<i<<endl;
//for (auto i : nodes2d) cout<<"nodi "<<i<<endl;
//for (auto i : weights2d) cout<<"pesi "<<i<<endl;
cout<<"The polygon is divided in these triangles:"<<endl;
auto V=Q.divide();
auto f=[](double x){return 0.0*x;};
unsigned int n=3; //number of points on each edge
cout<<"The global integral is = "<<Q.global_int(f,n)<<endl;

return 0;
}
