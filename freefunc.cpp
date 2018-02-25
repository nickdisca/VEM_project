#include "freefunc.hpp"

//computes the degrees in the order 00,10,01,20...
std::vector<std::array<int,2> > Polynomials(int k) {
	std::vector<std::array<int,2> > degree;
	for (int i=0; i<=k; i++) {
		for (int j=i; j>=0; j--){
			degree.push_back(std::array<int,2>{{j,i-j}});
		}
	}
return degree;
} 


//gives the GLL weights and points (naive way, can be made generic)
void computeDOF(std::vector<Point> const & points,unsigned int k,std::vector<double> & weights,std::vector<Point> & nodes){
	std::vector<double> W,N;
	double length; 
	nodes.clear(); weights.clear();
	if (k==1) {W.push_back(1.0); W.push_back(1.0);}
	if (k==2) {N.push_back(0.0); W.push_back(1.0/3.0); W.push_back(4.0/3.0); W.push_back(1.0/3.0);}
	if (k==3) {N.push_back(-1.0/std::sqrt(5.0)); N.push_back(1.0/std::sqrt(5.0)); 
		W.push_back(1.0/6.0); W.push_back(5.0/6.0); W.push_back(5.0/6.0); W.push_back(1.0/6.0);}
	if (k==4) {N.push_back(-std::sqrt(3.0/7.0)); N.push_back(0.0); N.push_back(std::sqrt(3.0/7.0)); 
		W.push_back(1.0/10.0); W.push_back(49.0/90.0); W.push_back(32.0/45.0); W.push_back(49.0/90.0); W.push_back(1.0/10.0);}
	if (k>=5) {std::cout<<"Error"<<std::endl; return;}

	//loop over all edges
	for (unsigned int i=0; i<points.size(); ++i){
		Point a=points[i], b=points[(i+1)%points.size()];
		length=distance(a,b);
		//map nodes and weights on each edge
		weights.push_back(length/2*W[0]);
		for (unsigned int j=0; j<N.size(); ++j) {
			nodes.push_back(Point((b[0]-a[0])/2*(1+N[j])+a[0],(b[1]-a[1])/2*(1+N[j])+a[1]));
			weights.push_back(length/2*W[j+1]);
		}
		weights.push_back(length/2*W[W.size()-1]); 
	}	

	return;
}
