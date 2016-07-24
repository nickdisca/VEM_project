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


return 0;
}
