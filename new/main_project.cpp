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

std::vector<double> nodes, weights;
double a=0.0,b=1.0;
unsigned int n=10;
Gauss_Leg(a,b,n,nodes,weights);
for (auto i : nodes) cout<<"nodi "<<i<<endl;
for (auto i : weights) cout<<"pesi "<<i<<endl;

Gauss_Leg(-1.0,4.0,6,nodes,weights);
for (auto i : nodes) cout<<"nodi "<<i<<endl;
for (auto i : weights) cout<<"pesi "<<i<<endl;

Gauss_Leg(1.0,2.0,4,nodes,weights);
for (auto i : nodes) cout<<"nodi "<<i<<endl;
for (auto i : weights) cout<<"pesi "<<i<<endl;

return 0;
}
