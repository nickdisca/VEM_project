#include <iostream>
#include <vector>
#include "Geometry.hpp"
#include "Mesh.hpp"

using namespace std;

int main()
{

double a=5.0,b=6.0;
Point c(a,b);
Point c2; c2.setCoordinates(b,a);
cout<<c<<c2<<c-c2<<c+c2<<3*c2;
cout<<"Distance "<<distance(c,c2)<<endl;

std::vector<Point> v;
v.push_back(Point(1.0,1.0)); v.push_back(Point{2.0,1.0}); v.push_back(Point(2.0,2.0)); v.push_back(Point(1.0,2.0));
//v.push_back(Point(0.0,0.0)); v.push_back(Point{1.0,0.0}); v.push_back(Point(1.0,1.0)); v.push_back(Point(0.0,1.0));

std::vector<unsigned int> uns; uns.push_back(0); uns.push_back(1); uns.push_back(2); uns.push_back(3);
Polygon p(uns,&v);
cout<<p;
cout<<"Area = "<<p.area()<<endl;
cout<<"Centroid is "<<p.centroid();
cout<<"Diameter = "<<p.diameter()<<endl;


return 0;
}
