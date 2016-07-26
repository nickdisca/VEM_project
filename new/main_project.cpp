#include <iostream>
#include <vector>
#include "Geometry.hpp"
#include "Mesh.hpp"

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
cout<<p;
cout<<p.ComputeD(2)<<endl;

return 0;
}
