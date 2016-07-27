#include "quadrature.hpp"
#include <cmath>
#include <limits>

std::ostream & operator << (std::ostream & ost, const Quadrature & Q){
	std::cout<<"Want to compute integrals on the following polygon:"<<std::endl<<Q.P<<std::endl;
	return ost;
}

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