#include "vem.hpp"

using namespace std;
using namespace Geometry;

Vertices AbstractPolygon::BoundaryDof(int k) const {
	Vertices vect(vertexes.size()*(k-1),0.0);
	if (k==1) return vect;
	double diffx{0},diffy{0};
	for (unsigned int i=1; i<vertexes.size()+1; i++) {
		diffx=((vertexes[i%vertexes.size()].x()-vertexes[i-1].x())/(k));
		diffy=((vertexes[i%vertexes.size()].y()-vertexes[i-1].y())/(k)); 
		
		for (int j=0; j<k-1; j++)
			vect[(i-1)*(k-1)+j].set(vertexes[i-1].x()+(j+1)*diffx,vertexes[i-1].y()+(j+1)*diffy);
	}
	return vect;
};

vector<array<int,2> > Polynomials(int k) {
	vector<array<int,2> > degree;
	for (int i=0; i<=k; i++) {
		for (int j=i; j>=0; j--){
			degree.push_back(array<int,2>{{j,i-j}});
			cout<<"{ "<<j<<" , "<<i-j<<" }"<<endl;
		}
	}
return degree;
} 

MyMatrix AbstractPolygon::ComputeD(int k) {
	MyMatrix D(vertexes.size()*k+k*(k-1)/2,(k+2)*(k+1)/2);
	cout<<"Created matrix D with size "<<D.GetRows()<<" x "<<D.GetCols()<<endl;

	Vertices BD=BoundaryDof(k);
	vector<array<int,2> > degree=Polynomials(k);
  	

	for (unsigned int j=0; j<D.GetCols(); j++) {
		for (unsigned int i=0; i<D.GetRows(); i++){
			array<int,2> actualdegree=degree[j];
			auto f=[=] (double x,double y)	
				{return pow((x-Centroid().x())/(Diameter()),actualdegree[0])*pow((y-Centroid().y())/(Diameter()),actualdegree[1]);};

			if (i<vertexes.size()) D(i,j)=f(vertexes[i].x(),vertexes[i].y());
			if (i>=vertexes.size() && i<vertexes.size()*2) {D(i,j)=f(BD[i-vertexes.size()].x(),BD[i-vertexes.size()].y());}
			if (i==vertexes.size()*2) D(i,j)=-1.0;
		}
	}


	return D;
}