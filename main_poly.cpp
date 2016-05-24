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

  //! five vertices
  Vertices v(5);

  v[0]={0.0,0.0};
  v[1]={3.0,0.0};
  v[2]={3.0,2.0};
  v[3]={1.5,4.0};
  v[4]={0.0,4.0};

  Polygon aPolygon(v);
  aPolygon.showMe();
  std::cout<<"Area: "<<aPolygon.area()<<std::endl;
  std::cout<<"Centroid: "<<aPolygon.Centroid().x()<<" "<<aPolygon.Centroid().y()<<std::endl;
  std::cout<<"Diameter: "<<aPolygon.Diameter()<<std::endl;

  /*
  //check that the default constructor works
  Grid F;
  cout<<"Numero di poligoni nella griglia: "<<F.grid_size()<<endl;
  cout<<"Area totale: "<<F.area()<<endl<<endl;

  //construction of a grid reading from the given file
  std::ifstream file("mesh.dat");
  Grid G(file);
  cout<<"Numero di poligoni nella griglia: "<<G.grid_size()<<endl;
  cout<<"Area totale: "<<G.area()<<endl<<endl;
  file.close();

  //check that the copy constructor works
  Grid H(G);
  cout<<"Numero di poligoni nella griglia: "<<H.grid_size()<<endl;
  cout<<"Area totale: "<<H.area()<<endl<<endl;

  //check that the assignment operator works
  F=H;
  cout<<"Numero di poligoni nella griglia: "<<F.grid_size()<<endl;
  cout<<"Area totale: "<<F.area()<<endl<<endl;

  //print the edges on the screen, just for simplicity
  F.printedges();

  //print the indexes of the edges in 3 separate files
  ofstream out1("AllEdges.dat");
  ofstream out2("BoundaryEdges.dat");
  ofstream out3("InternalEdges.dat");
  F.printedgesIndex(out1,out2,out3);
  */
}
  

