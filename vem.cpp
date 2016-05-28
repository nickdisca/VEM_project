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

double AbstractPolygon::IntegralWithDof(int k, int edge, int phi) {
	Vertices BD=BoundaryDof(k);
	Vertices nodes;
	nodes.push_back(vertexes[edge-1]);
	int j=0;
	for (unsigned int i=0; i<BD.size(); i++) 
		if ((edge-1)*(k-1)+j==i && j<=k-2) {nodes.push_back(BD[i]); j++;}
	nodes.push_back(vertexes[edge%vertexes.size()]);
	//for (auto i : nodes) cout<<"node: "<<i.x()<<" "<<i.y()<<endl;
	vector<double> weights(k+1,0.0);
	weights[0]=1.0/3/2; weights[1]=4.0/3/2; weights[2]=1.0/3/2;

	return weights[phi];
}

Point2D AbstractPolygon::Normal(int edge){
	Point2D N,aux=vertexes[edge]-vertexes[edge-1];
	N.set(aux.y(),-aux.x());
	//cout<<"The normal vector is: "<<N.x()<<" "<<N.y()<<endl;
	return N;
}

double AbstractPolygon::ComputeIntegral(int k, int d1, int d2) {
	double value{0.0};
	if (d1==0 && d2==0) return 1.0;
	if ((d1+d2)%2==1) return 0.0;
	if (d1==2 && d2==0) return 1.0/24;
	if (d1==0 && d2==2) return 1.0/24;
	if (d1==1 && d2==1) return 0.0;
	return value;
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
			if (i>=vertexes.size()*2) D(i,j)=ComputeIntegral(k,actualdegree[0],actualdegree[1]);
		}
	}


	return D;
}

MyMatrix AbstractPolygon::ComputeB(int k){
	MyMatrix B((k+2)*(k+1)/2,vertexes.size()*k+k*(k-1)/2);
	cout<<"Created matrix B with size "<<B.GetRows()<<" x "<<B.GetCols()<<endl;
	vector<array<int,2> > degree=Polynomials(k);
	B(0,B.GetCols()-1)=1.0;
	for (unsigned int j=0; j<B.GetCols()-1; j++) {
		for (unsigned int i=1; i<B.GetRows(); i++){
			array<int,2> actualdegree=degree[i];
			auto fx=[=] (double x,double y)	
{return actualdegree[0]/Diameter()*pow((x-Centroid().x())/(Diameter()),max(actualdegree[0]-1,0))*pow((y-Centroid().y())/(Diameter()),actualdegree[1]);};
			auto fy=[=] (double xx,double yy)	
{return actualdegree[1]/Diameter()*pow((xx-Centroid().x())/(Diameter()),actualdegree[0])*pow((yy-Centroid().y())/(Diameter()),max(actualdegree[1]-1,0));};
			

			B(i,j)=(Normal(j+1).x()*fx(vertexes[j].x(),vertexes[j].y())+
				Normal(j+1).y()*fy(vertexes[j].x(),vertexes[j].y()))*IntegralWithDof(k,j+1,0);
			cout<<"B("<<i<<","<<j<<")="<<B(i,j)<<endl;
			int next=(j==0 ? vertexes.size() : j);
			B(i,j)+=((Normal(next).x()*fx(vertexes[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y())+
				Normal(next).y()*fy(vertexes[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y()))*IntegralWithDof(k,next,2));
		
			cout<<"B("<<i<<","<<j<<")="<<B(i,j)<<endl;


		}
	}
	for (unsigned int i=0; i<B.GetRows(); i++) B(i,B.GetCols()-1)=-1.0;

return B;
}