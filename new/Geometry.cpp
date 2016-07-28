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

//calcola la normale al lato i con i che conta da 1 a N
Point Polygon::Normal(unsigned int edge_num){
	Point N;
	Point aux=(*pointer)[vertexes[edge_num%vertexes.size()]]-(*pointer)[vertexes[edge_num-1]];
	N.setCoordinates(aux[1],-aux[0]);
	//std::cout<<"The normal vector is: "<<N<<std::endl;
	N=N*(1.0/distance(Point(0.0,0.0),N));
	//std::cout<<"The normal versor is: "<<N<<std::endl;
	return N;
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
	double A(area()); std::cout<<"Area "<<A<<std::endl;
  	

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
			if (i>=vertexes.size()+dof.size()) {
				//std::cout<<i<<j<<std::endl;
				unsigned int ii=i-vertexes.size()-dof.size();
				Quadrature Q(*this);
				unsigned int pow1=actualdegree[0]+degree[ii][0], pow2=actualdegree[1]+degree[ii][1];
				//cambio la lambda perchè potrei avere int(m_alpha*polinomi)
				auto p= [C,diam,pow1,pow2] (double x,double y)
					{return pow((x-C[0])/diam,pow1)*pow((y-C[1])/diam,pow2);};

				//auto ppp=[=](double x,double y){return 1.0;}; D(i,j)=Q.global_int(ppp,k);
				D(i,j)=1.0/A*Q.global_int(p,k);
			}

		}
	}


	return D;
}



MatrixType Polygon::ComputeB(unsigned int k){

	//dimensions: nk x Ndof
	MatrixType B((k+2)*(k+1)/2,vertexes.size()*k+k*(k-1)/2);
	std::cout<<"Created matrix B with size "<<B.rows()<<" x "<<B.cols()<<std::endl;
	std::vector<Point> P=getPoints();
	//for (auto i: P) std::cout<<i;
	std::vector<Point> BD=getDof();
	//for (auto i: BD) std::cout<<i;
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter()); std::cout<<"Diameter "<<diam<<std::endl;
	Point C(centroid()); std::cout<<"Centroid "<<C;
	double A(area()); std::cout<<"Area "<<A<<std::endl;

	//ho bisogno di sapere i pesi associati (da migliorare assolutamente)
	std::vector<Point> dummy;
	std::vector<double> weights;
	std::cout<<"Here"<<std::endl;
	computeDOF(P,k,weights,dummy);
	//for (auto i : weights) std::cout<<i<<std::endl;


	unsigned int aux=0;
	for (unsigned int j=0; j<vertexes.size()+dof.size(); j++) {
		//da controllare posizione di aux
		int jj=(j-vertexes.size());
		if (j>=vertexes.size()){
		if (aux!=0 && aux%(k-1)==0) aux=aux+2;
		aux++;
		}

		for (unsigned int i=1; i<B.rows(); i++){
			std::array<int,2> actualdegree=degree[i];

			//calcolo le componenti del gradiente
			auto fx=[C,diam,actualdegree] (double x,double y)	
{return actualdegree[0]/diam*pow((x-C[0])/diam,std::max(actualdegree[0]-1,0))*pow((y-C[1])/diam,actualdegree[1]);};
			auto fy=[C,diam,actualdegree] (double xx,double yy)	
{return actualdegree[1]/diam*pow((xx-C[0])/diam,actualdegree[0])*pow((yy-C[1])/diam,std::max(actualdegree[1]-1,0));};
			

			std::cout<<i<<j<<std::endl;
			//contributi dovuti alle funzioni di base relative ai vertici
			if (j<vertexes.size()){
				B(i,j)=(Normal(j+1)[0]*fx(P[j][0],P[j][1])+Normal(j+1)[1]*fy(P[j][0],P[j][1]))*weights[j*(k+1)];
				//std::cout<<"B("<<i<<","<<j<<")="<<B(i,j)<<std::endl;
				unsigned int next=(j==0 ? vertexes.size() : j);
			//B(i,j)+=((Normal(next)[0]*fx(P[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y())+
			//	Normal(next).y()*fy(vertexes[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y()))*IntegralWithDof(k,next,2));
				unsigned int position=(j==0 ? weights.size()-1 : j*(k+1)-1);
				B(i,j)+=(Normal(next)[0]*fx(P[j][0],P[j][1])+Normal(next)[1]*fy(P[j][0],P[j][1]))*weights[position];
			}

			//contributi dovuti alle funzioni di base relative ai dof sul bordo
			else {
				//std::cout<<"Taking position number "<<aux<<std::endl;
				B(i,j)=(Normal(jj/(k-1)+1)[0]*fx(BD[jj][0],BD[jj][1])+
					Normal(jj/(k-1)+1)[1]*fy(BD[jj][0],BD[jj][1]))*weights[aux];
			}


		}
	}

//ultime colonne
for (unsigned int j=vertexes.size()+dof.size(); j<B.cols(); j++){
	for (unsigned int i=1; i<B.rows(); i++){
		std::array<int,2> actualdegree=degree[i];
		unsigned int jj=j-vertexes.size()-dof.size();
		std::cout<<i<<j<<std::endl;
		if (actualdegree[0]<=1 && actualdegree[1]<=1) B(i,j)=0.0;
		else {
		double coeff1=actualdegree[0]*(actualdegree[0]-1)/(diam*diam);
		double coeff2=actualdegree[1]*(actualdegree[1]-1)/(diam*diam);
		if (actualdegree[0]-2==degree[jj][0] && actualdegree[1]==degree[jj][1]) B(i,j)=-coeff1*A;
		if (actualdegree[1]-2==degree[jj][1] && actualdegree[0]==degree[jj][0]) B(i,j)=-coeff2*A;
		}
	}
}
//prima riga
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




//compute directly the matrix G (only for a control)
MatrixType Polygon::ComputeG(unsigned int k){

	//dimensions: nk x nk
	MatrixType G((k+2)*(k+1)/2,(k+2)*(k+1)/2);
	std::cout<<"Created matrix G with size "<<G.rows()<<" x "<<G.cols()<<std::endl;
	G.fill(0.0); //to be sure that it is initialized
	std::vector<Point> P=getPoints();
	//for (auto i: P) std::cout<<i;
	std::vector<Point> BD=getDof();
	//for (auto i: BD) std::cout<<i;
	std::vector<std::array<int,2> > degree=Polynomials(k);
	double diam(diameter()); std::cout<<"Diameter "<<diam<<std::endl;
	Point C(centroid()); std::cout<<"Centroid "<<C;
	double A(area()); std::cout<<"Area "<<A<<std::endl;

	//ho bisogno di sapere i pesi associati (da migliorare assolutamente)
	std::vector<Point> dummy;
	std::vector<double> weights;
	std::cout<<"Here"<<std::endl;
	computeDOF(P,k,weights,dummy);
	//for (auto i : weights) std::cout<<i<<std::endl;


	//ciclo sulle colonne (se parto da 0 calcola tutto lui, se parto da 1 sfrutto definizione)
	for (unsigned int j=1; j<G.cols(); j++) {


		for (unsigned int i=1; i<G.rows(); i++){
			std::array<int,2> actualdegreeI=degree[i];
			std::array<int,2> actualdegreeJ=degree[j];

			//calcolo le componenti del gradiente
		auto poli=[C,diam,actualdegreeJ] (double x,double y)	
{return pow((x-C[0])/diam,actualdegreeJ[0])*pow((y-C[1])/diam,actualdegreeJ[1]);};
			auto gradx=[C,diam,actualdegreeI] (double x,double y)	
{return actualdegreeI[0]/diam*pow((x-C[0])/diam,std::max(actualdegreeI[0]-1,0))*pow((y-C[1])/diam,actualdegreeI[1]);};
			auto grady=[C,diam,actualdegreeI] (double xx,double yy)	
{return actualdegreeI[1]/diam*pow((xx-C[0])/diam,actualdegreeI[0])*pow((yy-C[1])/diam,std::max(actualdegreeI[1]-1,0));};
			

			std::cout<<i<<j<<std::endl;

			for (unsigned int z=0; z<this->size(); z++) {
				//vertice iniziale di ogni lato
				//std::cout<<"vertex number"<<z<<"weight number"<<z*(k+1)<<std::endl;
				G(i,j)+=(Normal(z+1)[0]*gradx(P[z][0],P[z][1])+Normal(z+1)[1]*grady(P[z][0],P[z][1]))*
					poli(P[z][0],P[z][1])*weights[z*(k+1)];
				//vertice finale di ogni lato
				unsigned int aux=(z+1)%this->size();
				unsigned int position=(aux==0 ? weights.size()-1 : (k+1)*aux-1);
				//std::cout<<"vertex number"<<aux<<"weight number"<<position<<std::endl;
				G(i,j)+=(Normal(z+1)[0]*gradx(P[aux][0],P[aux][1])+Normal(z+1)[1]*grady(P[aux][0],P[aux][1]))*
					poli(P[aux][0],P[aux][1])*weights[position];
			}
			//BD intermedi
			int aux=-2;
			for (unsigned int z=0; z<BD.size(); z++){
				if (z%(k-1)==0) aux=aux+2;
				aux++;
				//std::cout<<"BD number"<<z<<"weight number"<<aux<<std::endl;
				G(i,j)+=(Normal(z/(k-1)+1)[0]*gradx(BD[z][0],BD[z][1])+
					Normal(z/(k-1)+1)[1]*grady(BD[z][0],BD[z][1]))*poli(BD[z][0],BD[z][1])*weights[aux];
			}

			//std::cout<<"i  "<<i<<" j   "<<j<<"    "<<G(i,j)<<std::endl;

			//termini interni
			if (actualdegreeI[0]<=1 && actualdegreeI[1]<=1) ;
			else {
			Quadrature Q(*this);
			//labmbda fun defining the integral function
			auto fun=[C,diam,actualdegreeI,actualdegreeJ] (double x, double y){
				int d1=actualdegreeI[0], d2=actualdegreeI[1];
				double res=d1*(d1-1)/(diam*diam)*pow((x-C[0])/diam,std::max(d1-2,0))*pow((y-C[1])/diam,d2);
				res+=d2*(d2-1)/(diam*diam)*pow((x-C[0])/diam,d1)*pow((y-C[1])/diam,std::max(d2-2,0));
				res=res*pow((x-C[0])/diam,actualdegreeJ[0])*pow((y-C[1])/diam,actualdegreeJ[1]);
				return res;
			}; 
			G(i,j)=G(i,j)-Q.global_int(fun,k);
			}

		} //ciclo su i
	} //ciclo su j


//prima riga
if (k>=2){
	for (unsigned int j=0; j<G.cols(); j++){
		std::array<int,2> grado=degree[j];
		auto poli=[C,diam,grado] (double x,double y)	
			{return pow((x-C[0])/diam,grado[0])*pow((y-C[1])/diam,grado[1]);};
		Quadrature Q(*this);
		G(0,j)=1.0/A*Q.global_int(poli,k);
	}
}
else {
	for (unsigned int j=0; j<G.cols(); j++){
		std::array<int,2> grado=degree[j];
		auto poli=[C,diam,grado] (double x,double y)	
			{return pow((x-C[0])/diam,grado[0])*pow((y-C[1])/diam,grado[1]);};
			for (unsigned int i=0; i<this->size(); i++) G(0,j)+=poli(P[i][0],P[i][1]);
		G(0,j)=G(0,j)/this->size();
	}
}


return G;
}