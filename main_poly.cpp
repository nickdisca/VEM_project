#include "Polygon.hpp"
#include "grid.hpp"
#include "edge.hpp"
#include <iostream>
#include <fstream>
//! Main program
int main()
{
  using namespace Geometry;
  using namespace std;
/*
  //! five vertices
  Vertices v(5);

  v[0]={0.5,0.5}; //C++11 sintax!
  v[1]={0.5,0.8};
  v[2]={0.0,0.0};
  v[3]={0.0,-1.0};
  v[4]={0.3,-0.5};

  Polygon aPolygon(v);
  aPolygon.showMe();
  std::cout<<"Area: "<<aPolygon.area()<<std::endl;
  // A triangle built with the first 3 vertices
  Triangle aTriangle(Vertices(v.begin(),v.begin()+3));
  AbstractPolygon * p_ab=&aTriangle;

  p_ab->showMe();
  std::cout<<"Area: "<<aTriangle.area()<<std::endl;
  //! Unit Square
  //C++11 syntax. I expoit the fact that the is an implicit conversion
  //between initializer lists with two double and a Point2d, since the latter has a  constructor taking two doubles
  //as argument 
  // In C++98 I would have written
  //Square aSquare(Point2D(0.0,0.0),1.0); 
  Square aSquare({0.0,0.0},1.0);
  Square s2(aSquare); 
  AbstractPolygon & r_ab=s2;
  r_ab.showMe();
  std::cout<<"Area: "<<r_ab.area()<<std::endl;
*/

  //construction of a grid reading from the given file
  std::ifstream file("mesh.dat");
  Grid G(file);
  cout<<"Numero di poligoni nella griglia: "<<G.grid_size()<<endl;
  cout<<"Area totale: "<<G.area()<<endl<<endl;
  file.close();

  //check that the copy constructor works
  //Grid H(G);
  //cout<<"Numero di poligoni nella griglia: "<<H.grid_size()<<endl;
  //cout<<"Area totale: "<<H.area()<<endl<<endl;

  //print the edges on the screen, just for simplicity
  G.printedges();

  //print the indexes of the edges in 3 separate files
  ofstream out1("AllEdges.dat");
  ofstream out2("BoundaryEdges.dat");
  ofstream out3("InternalEdges.dat");
  G.printedgesIndex(out1,out2,out3);

  return 0;
}
  

