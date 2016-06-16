#include "Polygon.hpp"
#include "grid.hpp"
#include "edge.hpp"
#include "vem.hpp"
#include "matrix.hpp"
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

  Vertices vv(4);

  vv[0]={0.0,0.0};
  vv[1]={1.0,0.0};
  vv[2]={1.0,1.0};
  vv[3]={0.0,1.0};

  Square aSquare(vv);
  aSquare.showMe();
  std::cout<<"Area: "<<aSquare.area()<<std::endl;
  std::cout<<"Centroid: "<<aSquare.Centroid().x()<<" "<<aSquare.Centroid().y()<<std::endl;
  std::cout<<"Diameter: "<<aSquare.Diameter()<<std::endl;

  MyMatrix A(5,5);
  A(1,1)=3.0; A(2,2)=5.0; A(3,3)=6.0;
  cout<<"My matrix is"<<endl<<A;

  std::vector<Point2D> BD1=aPolygon.BoundaryDof(1);
  std::cout<<"Case k=1:"<<std::endl;
  for (unsigned int i=0; i<BD1.size(); i++)
  	std::cout<<BD1[i].x()<<" "<<BD1[i].y()<<std::endl;

  std::vector<Point2D> BD2=aPolygon.BoundaryDof(2);
  std::cout<<"Case k=2:"<<std::endl;
  for (unsigned int i=0; i<BD2.size(); i++)
  	std::cout<<BD2[i].x()<<" "<<BD2[i].y()<<std::endl;

  std::vector<Point2D> BD3=aPolygon.BoundaryDof(3);
  std::cout<<"Case k=3:"<<std::endl;
  for (unsigned int i=0; i<BD3.size(); i++)
  	std::cout<<BD3[i].x()<<" "<<BD3[i].y()<<std::endl;

  std::vector<std::array<int,2> > degree1=Polynomials(1);
  std::cout<<std::endl;
  std::vector<std::array<int,2> > degree2=Polynomials(2);
  std::cout<<std::endl;
  std::vector<std::array<int,2> > degree3=Polynomials(3);
  std::cout<<std::endl;

  std::vector<Point2D> BDS2=aSquare.BoundaryDof(2);
  std::cout<<"Case k=2:"<<std::endl;
  for (unsigned int i=0; i<BDS2.size(); i++)
  std::cout<<BDS2[i].x()<<" "<<BDS2[i].y()<<std::endl;

  MyMatrix D=aSquare.ComputeD(2);
  std::cout<<"Matrix D: "<<std::endl<<D<<std::endl;

  MyMatrix B=aSquare.ComputeB(2);
  std::cout<<"Matrix B: "<<std::endl<<B<<std::endl;


  //construction of a grid reading from the given file
  std::ifstream file("squares.dat");
  Grid G(file);
  cout<<"Numero di poligoni nella griglia: "<<G.grid_size()<<endl;
  cout<<"Numero di elementi di bordo: "<<G.boundary_size()<<endl;
  cout<<"Numero di vertici totali: "<<G.vertices_size()<<endl;
  cout<<"Area totale: "<<G.area()<<endl<<endl;
  file.close();

  /*
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
  

