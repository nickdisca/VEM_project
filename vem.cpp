#include "vem.hpp"

using namespace std;

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