#include "quadrature.hpp"
#include <cmath>
#include <limits>

std::ostream & operator << (std::ostream & ost, const Quadrature & Q){
	std::cout<<"Want to compute integrals on the following polygon:"<<std::endl<<Q.P<<std::endl;
	return ost;
}

//calcola nodi e pesi di GL su un generico intervallo
void Gauss_Leg(double a, double b, unsigned int n, std::vector<double> & nodes,std::vector<double> & weights){
	double m=(n+1)/2.0;
	double xm=(a+b)/2,xl=(b-a)/2;
	double pi=4.0*std::atan(1.0);
	double pp=0.0,z=0.0;
	nodes.resize(n); weights.resize(n);

	for (unsigned int i=1; i<=m; i++){
		z=std::cos(pi*(i-0.25)/(n+0.5)); //std::cout<<z<<std::endl;
		//int cont=0;
		while(1){
			//cont++;
			double p1=1.0,p2=0.0;
			for (unsigned int j=1; j<=n; j++){
				double p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
				//std::cout<<p1<<p2<<p3<<std::endl;
			}
		pp=n*(z*p1-p2)/(z*z-1);
		//std::cout<<pp<<std::endl;
		double z1=z;
		z=z1-p1/pp;
		//std::cout<<"Hello"<<z<<"   "<<z1<<std::endl;
		if (std::abs(z-z1)<1000.0*std::numeric_limits<double>::epsilon()) {break;}
		}
		//std::cout<<"Here"<<i<<std::endl;
		nodes[i-1]=xm-z*xl;
		nodes[n-i]=xm+z*xl;
		weights[i-1]=2*xl/((1-z*z)*pp*pp);
		weights[n-i]=weights[i-1];
	}
	return;
}

//calcola nodi e pesi di GL 1D e 2D sul triangolo di riferimenti
//quelli 1D vanno bene solo sui cateti, ma potrebbero non essere necessari
void reference(unsigned int n,
	std::vector<double> & nodes1d, std::vector<double> & weights1d, std::vector<Point> & nodes2d, std::vector<double> & weights2d){
Gauss_Leg(0.0,1.0,n,nodes1d,weights1d);

std::vector<double> x,w;
Gauss_Leg(-1.0,1.0,n,x,w);

for (unsigned int i=0; i<n; i++){
	for (unsigned int j=0; j<n; j++){
		nodes2d.push_back(Point((1+x[i])/2.0,(1-x[i])*(1+x[j])/4.0));
		weights2d.push_back((1-x[i])*w[i]*w[j]/8.0);
	}
}
return;
}


std::vector<std::array<Point,3> > Quadrature::divide(){
	Point C=P.centroid();
	std::vector<Point> poi=P.getPoints();
	std::vector<std::array<Point,3> > vect;
	for (unsigned int i=0; i<poi.size(); i++) 
		vect.push_back({{poi[i],poi[(i+1)%P.size()],C}});
	for (auto i : vect) {for (auto j : i) std::cout<<j; std::cout<<std::endl;}
	return vect;
}

double Quadrature::local_int(std::array<Point,3> & tria, std::function<double(double,double)> f, unsigned int n){
	//std::cout<<"Hey"<<f(0.0)<<std::endl;
	std::vector<double> nodes1d, weights1d,weights2d;
	std::vector<Point> nodes2d;
	reference(n,nodes1d,weights1d,nodes2d,weights2d); //calcola nodi e pesi sul triangolo di riferimento

	MatrixType B,b;
	double det=map(tria,B,b);
	double res{0.0};

	for (unsigned int i=0; i<nodes2d.size(); i++){
		Point p=map_eval(nodes2d[i],B,b);
		res+=f(p[0],p[1])*weights2d[i];
	}
	//std::cout<<"Local integral"<<res*det<<std::endl;
	return res*det;
}

double Quadrature::global_int(std::function<double(double,double)> f, unsigned int n){
	double res{0.0};
	std::vector<std::array<Point,3> > triangles=this->divide();
	for (unsigned int i=0; i<triangles.size(); i++){
		res+=this->local_int(triangles[i],f,n);
	}
	return res;
}

double Quadrature::map(std::array<Point,3> & p, MatrixType & B, MatrixType & b){
	//MatrixType B(2,2, MatrixType b;
	B.resize(2,2); b.resize(2,1);
	Point aa=p[0],bb=p[1],cc=p[2];
	B(0,0)=bb[0]-aa[0]; B(0,1)=cc[0]-aa[0]; B(1,0)=bb[1]-aa[1]; B(1,1)=cc[1]-aa[1];
	b(0,0)=aa[0]; b(1,0)=aa[1];
	double det=B(1,1)*B(0,0)-B(1,0)*B(0,1);
	if (det==0) std::cout<<"Error: determinante = 0"<<std::endl;
	return std::abs(det);
}

Point Quadrature::map_eval(Point & x,MatrixType & B, MatrixType & b){
	MatrixType xx(2,1); xx(0,0)=x[0]; xx(1,0)=x[1];
	//MatrixType tmp=(B.lu()).solve(xx-b);
	MatrixType tmp=B*xx+b;
	return Point(tmp(0,0),tmp(1,0));
}