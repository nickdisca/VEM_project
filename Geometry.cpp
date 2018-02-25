#include <cmath>
#include <iostream>
#include "Geometry.hpp"
#include "quadrature.hpp"
#include "freefunc.hpp"

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
	ost<<"with area equal to: "<<p.area()<<std::endl<<"and centroid position: "<<p.centroid()<<std::endl;
	ost<<"with the following dofs: "<<std::endl;
	if (p.dof.size()==0) std::cout<<"Boundary dofs not set"<<std::endl;
	else
	for (unsigned int i=0; i<p.dof.size(); ++i) 
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
		Point const & p0(ver[(i-1+size) % size]);
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


//compute the normal versor w.r.t. edge i (i=1:n)
Point Polygon::Normal(unsigned int edge_num){

	Point N;

	//vector defining the edge
	Point aux=(*pointer)[vertexes[edge_num%vertexes.size()]]-(*pointer)[vertexes[edge_num-1]];

	//orthogonal vector
	N.setCoordinates(aux[1],-aux[0]);
	
	//normalize
	N=N*(1.0/distance(Point(0.0,0.0),N));
	
	return N;
}



MatrixType Polygon::ComputeD(unsigned int k) {
	
	//dimensions: Ndof x nk where Ndof=dimVk=nvert+nvert*(k-1)+n_(k-2)
	MatrixType D(vertexes.size()*k+k*(k-1)/2,(k+2)*(k+1)/2);
	
	//std::cout<<"Created matrix D with size "<<D.rows()<<" x "<<D.cols()<<std::endl; 
	
	D.fill(0.0);

	std::vector<Point> P=getPoints();
	std::vector<Point> BD=getDof();

	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());
	double A(area());
  	
	for (unsigned int j=0; j<D.cols(); j++) {
	
		std::array<int,2> actualdegree=degree[j];
		for (unsigned int i=0; i<D.rows(); i++){

			//polynomial to be evaluated
			auto f= [C,diam,actualdegree] (double x,double y)
				{return pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,actualdegree[1]);};

			//vertex
			if (i<vertexes.size()) 
				D(i,j)=f(P[i][0],P[i][1]);
				
			//GL point
			if (i>=vertexes.size() && i<vertexes.size()+dof.size()) {
				D(i,j)=f(BD[i-vertexes.size()][0],BD[i-vertexes.size()][1]);
			}
			
			//internal (integral of a polynomial)
			if (i>=vertexes.size()+dof.size()) {
				
				unsigned int ii=i-vertexes.size()-dof.size();
				
				Quadrature Q(*this);
				unsigned int pow1=actualdegree[0]+degree[ii][0], pow2=actualdegree[1]+degree[ii][1];
				
				//the polynomial to be integrated is m_alpha times a suitable monomial depending on dof
				auto p= [C,diam,pow1,pow2] (double x,double y)
					{return pow((x-C[0])/diam,pow1)*pow((y-C[1])/diam,pow2);};

				//compute integral using k 1D elements (enough to perform exact integration)
				D(i,j)=1.0/A*Q.global_int(p,k);
			}

		} //end loop i
	} //end loop j

	return D;
}



MatrixType Polygon::ComputeB(unsigned int k){

	//dimensions: nk x Ndof
	MatrixType B((k+2)*(k+1)/2,vertexes.size()*k+k*(k-1)/2);
	
	//std::cout<<"Created matrix B with size "<<B.rows()<<" x "<<B.cols()<<std::endl; 
	
	B.fill(0.0);
	std::vector<Point> P=getPoints();
	std::vector<Point> BD=getDof();
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());
	double A(area());

	//I need to know nodes and weights
	std::vector<Point> dummy;
	std::vector<double> weights;
	computeDOF(P,k,weights,dummy);
	
	int aux=0;
	
	for (unsigned int j=0; j<vertexes.size()+dof.size(); j++) {

		int jj=(j-vertexes.size());
		if (j>=vertexes.size()){
			if (aux!=0 && (aux+2)%(k+1)==0) //I need to skip weights if I find a vertex
				aux=aux+2;
			aux++;
		}

		for (unsigned int i=1; i<B.rows(); i++){
			
			std::array<int,2> actualdegree=degree[i];

			//components of the gradient
			auto fx=[C,diam,actualdegree] (double x,double y)	
				{return actualdegree[0]/diam*pow((x-C[0])/diam,std::max(actualdegree[0]-1,0))*pow((y-C[1])/diam,actualdegree[1]);};
			auto fy=[C,diam,actualdegree] (double xx,double yy)	
				{return actualdegree[1]/diam*pow((xx-C[0])/diam,actualdegree[0])*pow((yy-C[1])/diam,std::max(actualdegree[1]-1,0));};
			
			//vertex: contribution from two edges
			if (j<vertexes.size()){
				B(i,j)=(Normal(j+1)[0]*fx(P[j][0],P[j][1])+Normal(j+1)[1]*fy(P[j][0],P[j][1]))*weights[j*(k+1)];
				
				unsigned int next=(j==0 ? vertexes.size() : j); //which is the second edge
				unsigned int position=(j==0 ? weights.size()-1 : j*(k+1)-1); //next weight position
			
				B(i,j)+=(Normal(next)[0]*fx(P[j][0],P[j][1])+Normal(next)[1]*fy(P[j][0],P[j][1]))*weights[position];
			}

			//GL nodes
			else {
				B(i,j)=(Normal(jj/(k-1)+1)[0]*fx(BD[jj][0],BD[jj][1])+
					Normal(jj/(k-1)+1)[1]*fy(BD[jj][0],BD[jj][1]))*weights[aux];
			}


		} //end loop i
	} //end loop j


	//last columns (internal dofs)
	for (unsigned int j=vertexes.size()+dof.size(); j<B.cols(); j++){
		for (unsigned int i=1; i<B.rows(); i++){

			std::array<int,2> actualdegree=degree[i];
			unsigned int jj=j-vertexes.size()-dof.size();

			if (actualdegree[0]<=1 && actualdegree[1]<=1) ;//do nothing (Laplace op is zero)
			else {
				double coeff1=actualdegree[0]*(actualdegree[0]-1)/(diam*diam);
				double coeff2=actualdegree[1]*(actualdegree[1]-1)/(diam*diam);
				if (actualdegree[0]-2==degree[jj][0] && actualdegree[1]==degree[jj][1]) 
					B(i,j)+=-coeff1*A;
				else 
					if (actualdegree[1]-2==degree[jj][1] && actualdegree[0]==degree[jj][0]) 
						B(i,j)+=-coeff2*A;
				}
		} //end loop i
	} //end loop j


	//first row (different definition of k=1 and k>=2)
	if (k>=2){
		for (unsigned int j=0; j<B.cols(); j++)
			B(0,j)=(j==vertexes.size()+dof.size() ? 1.0 : 0.0);
	}
	else {
		for (unsigned int j=0; j<B.cols(); j++)
			B(0,j)=(j<vertexes.size()+dof.size() ? 1.0/vertexes.size() : 0.0);
	}

	return B;
}




//compute directly the matrix G (not necessarily required, but it is a useful check)
//similar strategy w.r.t. matrix B
//alternatively: integrate directly the inner product
MatrixType Polygon::ComputeG(unsigned int k){

	//dimensions: nk x nk
	MatrixType G((k+2)*(k+1)/2,(k+2)*(k+1)/2);
	
	//std::cout<<"Created matrix G with size "<<G.rows()<<" x "<<G.cols()<<std::endl; 
	
	G.fill(0.0);
	std::vector<Point> P=getPoints();
	std::vector<Point> BD=getDof();
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());
	double A(area());

	std::vector<Point> dummy;
	std::vector<double> weights;
	computeDOF(P,k,weights,dummy);

	//note that the first column is always zero, because gradient(constant)=0
	for (unsigned int j=1; j<G.cols(); j++) {

		for (unsigned int i=1; i<G.rows(); i++){
			std::array<int,2> actualdegreeI=degree[i];
			std::array<int,2> actualdegreeJ=degree[j];

			//polynomial and components of the gradient
			auto poli=[C,diam,actualdegreeJ] (double x,double y)	
				{return pow((x-C[0])/diam,actualdegreeJ[0])*pow((y-C[1])/diam,actualdegreeJ[1]);};
			auto gradx=[C,diam,actualdegreeI] (double x,double y)	
				{return actualdegreeI[0]/diam*pow((x-C[0])/diam,std::max(actualdegreeI[0]-1,0))*pow((y-C[1])/diam,actualdegreeI[1]);};
			auto grady=[C,diam,actualdegreeI] (double xx,double yy)	
				{return actualdegreeI[1]/diam*pow((xx-C[0])/diam,actualdegreeI[0])*pow((yy-C[1])/diam,std::max(actualdegreeI[1]-1,0));};
			
			//compute boundary integral (contribution from all dofs)
			//vertexes
			for (unsigned int z=0; z<this->size(); z++) {
				
				//initial point
				G(i,j)+=(Normal(z+1)[0]*gradx(P[z][0],P[z][1])+Normal(z+1)[1]*grady(P[z][0],P[z][1]))*
					poli(P[z][0],P[z][1])*weights[z*(k+1)];
				
				//end point
				unsigned int aux=(z+1)%this->size();
				unsigned int position=(aux==0 ? weights.size()-1 : (k+1)*aux-1);

				G(i,j)+=(Normal(z+1)[0]*gradx(P[aux][0],P[aux][1])+Normal(z+1)[1]*grady(P[aux][0],P[aux][1]))*
					poli(P[aux][0],P[aux][1])*weights[position];
			}
			
			//GL points
			int aux=-2;
			for (unsigned int z=0; z<BD.size(); z++){
				if (z%(k-1)==0) 
					aux=aux+2;
				
				aux++;
			
				G(i,j)+=(Normal(z/(k-1)+1)[0]*gradx(BD[z][0],BD[z][1])+
					Normal(z/(k-1)+1)[1]*grady(BD[z][0],BD[z][1]))*poli(BD[z][0],BD[z][1])*weights[aux];
			}


			//compute internal integral
			if (actualdegreeI[0]<=1 && actualdegreeI[1]<=1) ; //do nothing
			else {
				Quadrature Q(*this);
				//labmbda function defining the integral function
				auto fun=[C,diam,actualdegreeI,actualdegreeJ] (double x, double y){
					int d1=actualdegreeI[0], d2=actualdegreeI[1];
					double res=d1*(d1-1)/(diam*diam)*pow((x-C[0])/diam,std::max(d1-2,0))*pow((y-C[1])/diam,d2);
					res+=d2*(d2-1)/(diam*diam)*pow((x-C[0])/diam,d1)*pow((y-C[1])/diam,std::max(d2-2,0));
					res=res*pow((x-C[0])/diam,actualdegreeJ[0])*pow((y-C[1])/diam,actualdegreeJ[1]);
					return res;
				};
			//compute integral with k 1D points (enough to perform exact integration) 
			G(i,j)=G(i,j)-Q.global_int(fun,k);
			}

		} //end loop i
	} //end loop j


	//first row
	for (unsigned int j=0; j<G.cols(); j++){
		std::array<int,2> deg=degree[j];
		auto poli=[C,diam,deg] (double x,double y)	
			{return pow((x-C[0])/diam,deg[0])*pow((y-C[1])/diam,deg[1]);};

		if (k>=2) {
			Quadrature Q(*this);
			G(0,j)=1.0/A*Q.global_int(poli,k);
		}
		else {
			for (unsigned int i=0; i<this->size(); i++) 
				G(0,j)+=poli(P[i][0],P[i][1]);
				
			G(0,j)=G(0,j)/this->size();
		}
	}

	return G;
}


//stiffness matrix for the Laplace problem
MatrixType Polygon::LocalStiffness(unsigned int k){
	MatrixType B=ComputeB(k), D=ComputeD(k);
	MatrixType G=ComputeG(k);
	MatrixType Pistar=(G.lu()).solve(B);
	MatrixType Pi=D*Pistar;
	MatrixType I; I.setIdentity(Pi.rows(),Pi.cols());

	//now we define G_tilde, which is obtained from G by setting the first row to zero
	G.topRows(1).fill(0.0);

	return Pistar.transpose()*G*Pistar+(I-Pi).transpose()*(I-Pi);
}

//stiffness (diffusion) matrix for a generic elliptic problem
MatrixType Polygon::LocalStiffness_weighted(unsigned int k, std::function<double (double,double)> mu, double mu_bar, bool constant_mu){

	//if the viscosity is constant, use standard algorithm
	if (constant_mu) {
		//std::cout<<"Computing constant stiffness matrix"<<std::endl;
		return mu(0,0)*LocalStiffness(k);
	}

	//std::cout<<"Computing weighted stiffness matrix"<<std::endl;
	
	MatrixType B=ComputeB(k), D=ComputeD(k);
	MatrixType G=ComputeG(k);
	MatrixType Pi=D*(G.lu()).solve(B);
	MatrixType I; I.setIdentity(Pi.rows(),Pi.cols());
	MatrixType H=ComputeH(k-1,[](double x, double y){return 1.0;});
	MatrixType H_weight=ComputeH(k-1,mu);
	MatrixType Ex=ComputeE(k,0);
	MatrixType Ey=ComputeE(k,1);
	
	MatrixType Pi0_starx=(H.lu()).solve(Ex);
	MatrixType Pi0_stary=(H.lu()).solve(Ey);
	
	return Pi0_starx.transpose()*H_weight*Pi0_starx+Pi0_stary.transpose()*H_weight*Pi0_stary+mu_bar*(I-Pi).transpose()*(I-Pi);

}




MatrixType Polygon::LocalTransport(unsigned int k, std::function<double (double,double)> beta_x,std::function<double (double,double)> beta_y) {
	
	MatrixType B=ComputeB(k), D=ComputeD(k);
	MatrixType G=ComputeG(k);

	MatrixType H_minus=ComputeH(k-1,[](double x, double y){return 1.0;});
	MatrixType H=ComputeH(k,[](double x, double y){return 1.0;});
	MatrixType H_x=ComputeH(k,beta_x,k-1,k);
	MatrixType H_y=ComputeH(k,beta_y,k-1,k);
	MatrixType Ex=ComputeE(k,0);
	MatrixType Ey=ComputeE(k,1);
	
	MatrixType Pi0_starx=(H_minus.lu()).solve(Ex);
	MatrixType Pi0_stary=(H_minus.lu()).solve(Ey);
	
	MatrixType C=ComputeC(k);
	MatrixType Pi0_star=(H.lu()).solve(C);

	MatrixType transport=Pi0_starx.transpose()*H_x*Pi0_star+Pi0_stary.transpose()*H_y*Pi0_star;

	return transport.transpose();
}




MatrixType Polygon::ComputeH(unsigned int k, std::function<double (double,double)> weight,
	unsigned int krows, unsigned int kcols){
	
	//dimensions: nkrows x nkcols (coincide except for transport term which have k-1 and k...)	
	MatrixType H((krows+1)*(krows+2)/2,(kcols+1)*(kcols+2)/2);
	
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());

	for (unsigned int i=0; i<H.rows(); i++){
		for (unsigned int j=0; j<H.cols(); j++){
			
			std::array<int,2> dgI=degree[i];
			std::array<int,2> dgJ=degree[j];
			Quadrature Q(*this);
			
			//compute the product of the polynomials and use the weight
			auto poli=[C,diam,dgI,dgJ,weight](double x,double y)
			{double res=pow((x-C[0])/diam,dgI[0])*pow((y-C[1])/diam,dgI[1]);
				res=res*pow((x-C[0])/diam,dgJ[0])*pow((y-C[1])/diam,dgJ[1]);
				res=res*weight(x,y);
				return res;};
				
			//compute the integral
			//note: the weight can be a generic function, I may need to overintegrate 
			//for our tests k+3 is enough, but it can be changed
			H(i,j)=Q.global_int(poli,k+3);
			
		} //end loop j
	} //end loop i
	
	return H;
}


MatrixType Polygon::ComputeC(unsigned int k){
	
	MatrixType C((k+1)*(k+2)/2,vertexes.size()*k+k*(k-1)/2);
	
	//std::cout<<"Created matrix C with size "<<C.rows()<<" x "<<C.cols()<<std::endl;
	
	std::vector<std::array<int,2> > degree=Polynomials(k);
	Point centr(centroid());
	double A(area());

	MatrixType M=ComputeH(k)*((ComputeG(k).lu()).solve(ComputeB(k)));
	
	for (unsigned int alpha=0; alpha<C.rows(); alpha++){
		for (unsigned int j=0; j<C.cols(); j++){
			
			int jj=j-vertexes.size()-dof.size();
			
			//I can compute the L^2 projection using the dofs
			if(alpha<k*(k-1)/2) { 
					if (jj==(int)alpha) 
						C(alpha,j)=A;
					else C(alpha,j)=0.0;
			}
			
			//I need to use the Pi^\nabla projection
			else 
				C(alpha,j)=M(alpha,j);
				
		} //end loop j
	} //end loop alpha
	
	return C;
}


MatrixType Polygon::LoadTerm(unsigned int k, std::function<double (double,double)> f){
	
	//dimension: ndof
	MatrixType F(vertexes.size()*k+k*(k-1)/2,1); 
	
	F.fill(0.0);
	//std::cout<<"Created vector F with size "<<F.rows()<<std::endl;
		
	MatrixType M=(ComputeH(k).lu()).solve(ComputeC(k));
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());

	for (unsigned int i=0; i<F.rows(); i++){
		
		//I need to integrate f*m_alpha
		for (unsigned int alpha=0; alpha<(k+1)*(k+2)/2; alpha++){
			
			Quadrature Q(*this);
			std::array<int,2> actualdegree=degree[alpha];
			
			auto fun= [actualdegree,C,diam,f](double x,double y) {
				return f(x,y)*pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,actualdegree[1]);};

			//compute the integrals: for our tests k points are enough, but it can be modified
			double integral=Q.global_int(fun,k);

			F(i,0)+=M(alpha,i)*integral;
		}

	} //end loop i
	return F;
}


MatrixType Polygon::LocalConvert(unsigned int k, std::function<double (double,double)> uex){
	
	//dimension: ndof
	MatrixType U(vertexes.size()*k+k*(k-1)/2,1); 
	
	U.fill(0.0);

	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter());
	Point C(centroid());
	double A(area());
	std::vector<Point> P=getPoints();
	std::vector<Point> BD=getDof();

	for (unsigned int i=0; i<U.rows(); i++){
		
		//vertexes (function evaluation)
		if (i<vertexes.size()) 
			U(i,0)=uex(P[i][0],P[i][1]);
		
		//GL points (function evaluation)
		if (i>=vertexes.size() && i<k*vertexes.size()) 
			U(i,0)=uex(BD[i-vertexes.size()][0],BD[i-vertexes.size()][1]);
		
		//internal dof (compute integral)
		if (i>=k*vertexes.size()) {
				unsigned int ii=i-k*vertexes.size();
				Quadrature Q(*this);
				std::array<int,2> actualdegree=degree[ii];
				auto fun= [actualdegree,C,diam,uex](double x,double y) {
					return uex(x,y)*pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,actualdegree[1]);};

			//I may need to change the number of quadrature points
			U(i,0)=1.0/A*Q.global_int(fun,k);
			}

	} //end loop i
	
	return U;
}


MatrixType Polygon::LocalMass(unsigned int k){
	MatrixType C=ComputeC(k), H=ComputeH(k), D=ComputeD(k);
	MatrixType Pi0=D*(H.lu()).solve(C);
	MatrixType I; I.setIdentity(Pi0.rows(),Pi0.cols());

	return C.transpose()*(H.lu()).solve(C)+area()*(I-Pi0).transpose()*(I-Pi0);
}


//same strategy as matrix B
//VAR denotes the variable w.r.t. I am differentiating (VAR=0 is d/dx, VAR=1 is d/dy)
MatrixType Polygon::ComputeE(unsigned int k, unsigned int VAR){

	//dimensions: nk x ndof
	MatrixType E((k+1)*(k)/2,vertexes.size()*k+k*(k-1)/2);

	//std::string output=(VAR==0)?"x":"y";
	//std::cout<<"Computing matrix associated with derivative with respect to "<<output<<std::endl;

	E.fill(0.0);
	std::vector<Point> P=getPoints();
	std::vector<Point> BD=getDof();
	std::vector<std::array<int,2> > degree=Polynomials(k-1);
	double diam(diameter());
	Point C(centroid());
	double A(area());

	std::vector<Point> dummy;
	std::vector<double> weights;
	computeDOF(P,k,weights,dummy);
	
	int aux=0;
	for (unsigned int j=0; j<vertexes.size()+dof.size(); j++) {
		
		int jj=(j-vertexes.size());
		if (j>=vertexes.size()){
			if (aux!=0 && (aux+2)%(k+1)==0) 
				aux=aux+2;
			aux++;
		}

		for (unsigned int i=0; i<E.rows(); i++){
			std::array<int,2> actualdegree=degree[i];

			//compute the polynomial and its gradient
			auto poli=[C,diam,actualdegree] (double x, double y)
				{ return pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,actualdegree[1]); };
			auto der_x=[C,diam,actualdegree] (double x,double y)	
				{return actualdegree[0]/diam*pow((x-C[0])/diam,std::max(actualdegree[0]-1,0))*pow((y-C[1])/diam,actualdegree[1]);};
			auto der_y=[C,diam,actualdegree] (double x,double y)
				{return actualdegree[1]/diam*pow((x-C[0])/diam,actualdegree[0])*pow((y-C[1])/diam,std::max(actualdegree[1]-1,0));};

	
			//vertexes
			if (j<vertexes.size()){
				E(i,j)=Normal(j+1)[VAR]*poli(P[j][0],P[j][1])*weights[j*(k+1)];
				
				unsigned int next=(j==0 ? vertexes.size() : j);
				unsigned int position=(j==0 ? weights.size()-1 : j*(k+1)-1);
				
				E(i,j)+=Normal(next)[VAR]*poli(P[j][0],P[j][1])*weights[position];
			}

			//GL points
			else 
				E(i,j)=Normal(jj/(k-1)+1)[VAR]*poli(BD[jj][0],BD[jj][1])*weights[aux];


		} //end loop i
	} //end loop j

	//last columns
	for (unsigned int j=vertexes.size()+dof.size(); j<E.cols(); j++){
		for (unsigned int i=0; i<E.rows(); i++){
			std::array<int,2> actualdegree=degree[i];
			unsigned int jj=j-vertexes.size()-dof.size();

			double coeff=actualdegree[VAR]/diam;
			if (VAR==0 && actualdegree[0]-1==degree[jj][0] && actualdegree[1]==degree[jj][1]) 
				E(i,j)+=-coeff*A;
			if (VAR==1 && actualdegree[0]==degree[jj][0] && actualdegree[1]-1==degree[jj][1])
				E(i,j)+=-coeff*A;
			
		}
	}

	return E;
}



