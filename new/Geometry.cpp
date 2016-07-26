#include <cmath>
#include <iostream>
#include "Geometry.hpp"

Point Point::operator *(const double & d) const
  {return Point(d*M_coor[0],d*M_coor[1]);};

double distance(const Point & a, const Point & b){
	Point tmp(b-a);
	return std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
};

std::ostream & operator << (std::ostream & ost, const Point & p){
	ost<<"X = "<<p.M_coor[0]<<" Y = "<<p.M_coor[1]<<std::endl;
	return ost;
};




std::vector<Point> Polygon::getPoints() const{
	std::vector<Point> tmp(this->size());
	for (unsigned int i=0; i<tmp.size(); ++i) 
		tmp[i]=pointer->operator[](vertexes[i]);
	return tmp;
};


std::ostream & operator << (std::ostream & ost, const Polygon & p){
	ost<<"The polygon is the following: "<<std::endl;
	for (unsigned int i=0; i<p.size(); ++i) 
		ost<<"Index "<<" = "<<p.vertexes[i]<<" corresponding to point "<<p.pointer->operator[](p.vertexes[i]);
	ost<<"with the following dofs: "<<std::endl;
	if (p.dof.size()==0) std::cout<<"Boundary dofs not set"<<std::endl;
	else
	for (unsigned int i=0; i<p.size(); ++i) 
		ost<<"Index "<<" = "<<p.dof[i]<<" corresponding to point "<<p.pointer_dof->operator[](p.dof[i]);
return ost;
};


double Polygon::area() const{
	auto size=this->size();
	if (size<3) return 0.0;
	double result(0);
	auto ver(this->getPoints());

	for (decltype(size) i=0; i<size;++i){
		Point const & p1(ver[i]);
		Point const & p2(ver[(i+1) % size]);
		Point const & p0(ver[(i-1) % size]);
		result+=p1[0]*(p2[1]-p0[1]);
	}
	return 0.5*result;
};

Point Polygon::centroid() const{
	double x=0.0,y=0.0;
	auto ver(this->getPoints());
	unsigned int size=ver.size();
	for (unsigned int i=0; i<size; ++i){
	  x+=(ver[i][0]+ver[(i+1)%size][0])*(ver[i][0]*ver[(i+1)%size][1]-ver[(i+1)%size][0]*ver[i][1]);
	  y+=(ver[i][1]+ver[(i+1)%size][1])*(ver[i][0]*ver[(i+1)%size][1]-ver[(i+1)%size][0]*ver[i][1]);
	}
	x=x/(6*area());
	y=y/(6*area());
	return Point{x,y};
};

double Polygon::diameter() const{
	double d{0.0};
	auto ver(this->getPoints());
	for (unsigned int i=0; i<this->size()-1; i++){
		for (unsigned int j=0; j<this->size(); j++)
			d=std::max(d,distance(ver[i],ver[j]));
	}
	return d;
};

void Polygon::setDof(std::vector<unsigned int> const & v, std::vector<Point> * p){
	dof=v; pointer_dof=p;
	return;
}

std::vector<Point> Polygon::getDof() const{
	std::vector<Point> tmp(dof.size());
	for (unsigned int i=0; i<tmp.size(); ++i) 
		tmp[i]=pointer_dof->operator[](dof[i]);
	return tmp;
};

//calcola i gradi dei polinomi nell'ordine 00,10,01,20...
std::vector<std::array<int,2> > Polynomials(int k) {
	std::vector<std::array<int,2> > degree;
	for (int i=0; i<=k; i++) {
		for (int j=i; j>=0; j--){
			degree.push_back(std::array<int,2>{{j,i-j}});
			//std::cout<<"{ "<<j<<" , "<<i-j<<" }"<<std::endl;
		}
	}
return degree;
} 

MatrixType Polygon::ComputeD(unsigned int k) {
	
	//dimensions: Ndof x nk dove Ndof=dimVk=nvert+nvert*(k-1)+n_(k-2)
	MatrixType D(vertexes.size()*k+k*(k-1)/2,(k+2)*(k+1)/2);
	std::cout<<"Created matrix D with size "<<D.rows()<<" x "<<D.cols()<<std::endl;

	std::vector<Point> P=getPoints();
	//for (auto i: P) std::cout<<i;
	std::vector<Point> BD=getDof();
	//for (auto i: BD) std::cout<<i;
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter()); std::cout<<"Diameter "<<diam<<std::endl;
	Point C(centroid()); std::cout<<"Centroid "<<C;
  	

	for (unsigned int j=0; j<D.cols(); j++) {
		std::array<int,2> actualdegree=degree[j];
		for (unsigned int i=0; i<D.rows(); i++){

			//polinomio da valutare
			auto f= [C,diam,actualdegree] (double x,double y)
				{return pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,actualdegree[1]);};

			//inserisco nella matrice
			if (i<vertexes.size()) D(i,j)=f(P[i][0],P[i][1]);
			if (i>=vertexes.size() && i<vertexes.size()+dof.size()) 
			{
				D(i,j)=f(BD[i-vertexes.size()][0],BD[i-vertexes.size()][1]);
				//std::cout<<i<<j<<"   "<<BD[i-vertexes.size()][0]<<"  "<<BD[i-vertexes.size()][1]<<std::endl;
			}
			//devo calcolare gli integrali dei polinomi
			//if (i>=vertexes.size()+dof.size()) D(i,j)=ComputeIntegral(k,actualdegree[0],actualdegree[1]);
		}
	}


	return D;
}