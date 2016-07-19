#include "vem.hpp"

using namespace std;
using namespace Geometry;
using MatrixType=Matrix<double,Dynamic,Dynamic>;

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
	if (d1==1 && d2==0) return 0.0;
	if (d1==0 && d2==1) return 0.0;
	if (d1==2 && d2==0) return 1.0/24;
	if (d1==0 && d2==2) return 1.0/24;
	if (d1==1 && d2==1) return 0.0;
	if (d1==3 && d2==0) return 0.0;
	if (d1==2 && d2==1) return 0.0;
	if (d1==1 && d2==2) return 0.0;
	if (d1==0 && d2==3) return 0.0;
	if (d1==4 && d2==0) return 1.0/320;
	if (d1==3 && d2==1) return 0.0;
	if (d1==2 && d2==2) return 1.0/576;
	if (d1==1 && d2==3) return 0.0;
	if (d1==0 && d2==4) return 1.0/320;
	cout<<"Unknown value"<<endl;
	return value;
}

MatrixType AbstractPolygon::ComputeD(int k) {
	MatrixType D(vertexes.size()*k+k*(k-1)/2,(k+2)*(k+1)/2);
	cout<<"Created matrix D with size "<<D.rows()<<" x "<<D.cols()<<endl;

	Vertices BD=BoundaryDof(k);
	vector<array<int,2> > degree=Polynomials(k);
  	

	for (unsigned int j=0; j<D.cols(); j++) {
		for (unsigned int i=0; i<D.rows(); i++){
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

MatrixType AbstractPolygon::ComputeH(int k){
MatrixType H((k+2)*(k+1)/2,(k+2)*(k+1)/2);
vector<array<int,2> > degree=Polynomials(k);
//devo calcolare gli integrali dei prodotti dei monomi di base
for (unsigned int i=0; i<H.rows(); i++){
	for (unsigned int j=0; j<H.cols(); j++){
		H(i,j)=ComputeIntegral(k,degree[i][0]+degree[j][0],degree[i][1]+degree[j][1]);
	}
}
return H;
}

MatrixType AbstractPolygon::LoadTerm(int k){
	MatrixType F(vertexes.size()*k+k*(k-1)/2,1); //suppose load term constant=1
	vector<array<int,2> > degree=Polynomials(k);
	MatrixType Pi0star=((ComputeH(k).lu()).solve(ComputeC(k)));
	cout<<Pi0star<<endl;
	for (unsigned int i=0; i<F.rows(); i++){
		for (unsigned int alpha=0; alpha<(k+2)*(k+1)/2; alpha++){
			F(i,0)+=Pi0star(alpha,i)*ComputeIntegral(k,degree[alpha][0],degree[alpha][1]);
			//cout<<F(i,0)<<"  "<<i<<endl;
		}
	}
	return F;
}

MatrixType AbstractPolygon::ComputeC(int k){
MatrixType C((k+1)*(k+2)/2,vertexes.size()*k+k*(k-1)/2);
for (unsigned int alpha=0; alpha<C.rows(); alpha++){
	for (unsigned int i=0; i<C.cols(); i++){
		if (alpha<k*(k-1)/2 && i==k*vertexes.size()+alpha) C(alpha,i)=area();
		if (alpha<k*(k-1)/2 && i!=k*vertexes.size()+alpha) C(alpha,i)=0.0;
		if (alpha>= k*(k-1)/2) {MatrixType M= ComputeH(k)*((ComputeG(k).lu()).solve(ComputeB(k)));
							C(alpha,i)=M(alpha,i);}
	}
}
return C;
}

MatrixType AbstractPolygon::ComputeB(int k){
	MatrixType B((k+2)*(k+1)/2,vertexes.size()*k+k*(k-1)/2);
	cout<<"Created matrix B with size "<<B.rows()<<" x "<<B.cols()<<endl;
	Vertices BD=BoundaryDof(k);
	vector<array<int,2> > degree=Polynomials(k);

	for (unsigned int j=0; j<B.cols()-1; j++) {
		for (unsigned int i=1; i<B.rows(); i++){
			array<int,2> actualdegree=degree[i];
			auto fx=[=] (double x,double y)	
{return actualdegree[0]/Diameter()*pow((x-Centroid().x())/(Diameter()),max(actualdegree[0]-1,0))*pow((y-Centroid().y())/(Diameter()),actualdegree[1]);};
			auto fy=[=] (double xx,double yy)	
{return actualdegree[1]/Diameter()*pow((xx-Centroid().x())/(Diameter()),actualdegree[0])*pow((yy-Centroid().y())/(Diameter()),max(actualdegree[1]-1,0));};
			
			if (j<vertexes.size()){
			B(i,j)=(Normal(j+1).x()*fx(vertexes[j].x(),vertexes[j].y())+Normal(j+1).y()*fy(vertexes[j].x(),vertexes[j].y()))*IntegralWithDof(k,j+1,0);
			//cout<<"B("<<i<<","<<j<<")="<<B(i,j)<<endl;
			int next=(j==0 ? vertexes.size() : j);
			B(i,j)+=((Normal(next).x()*fx(vertexes[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y())+
				Normal(next).y()*fy(vertexes[next%(vertexes.size())].x(),vertexes[next%(vertexes.size())].y()))*IntegralWithDof(k,next,2));
			}
			else {
				double jj=j%vertexes.size();
			B(i,j)=(Normal(jj+1).x()*fx(BD[jj].x(),BD[jj].y())+Normal(jj+1).y()*fy(BD[jj].x(),BD[jj].y()))*
			IntegralWithDof(k,jj+1,1);
			}


		}
	}

	B(0,B.cols()-1)=1.0;
	for (unsigned int i=1; i<B.rows(); i++) {
		array<int,2> actualdegree=degree[i];
		if (actualdegree[0]<=1 || actualdegree[1]<=1) B(i,B.cols()-1)=0.0;
		if (i==3 || i==5) B(i,B.cols()-1)=-area()*2/2;
	}

return B;
}

MatrixType AbstractPolygon::ComputeG(int k){
	return ComputeB(k)*ComputeD(k);
}

MatrixType AbstractPolygon::ComputeGTilda(int k){
	MatrixType G=ComputeG(k);
	for (unsigned int j=0; j<G.cols(); j++) G(0,j)=0.0;
	return G;
}

MatrixType AbstractPolygon::ComputeStiffness(int k){
	MatrixType Pistar=(ComputeG(k).lu()).solve(ComputeB(k)); //cout<<"Pi star"<<Pistar.rows()<<Pistar.cols()<<endl;
	MatrixType Pi=ComputeD(k)*Pistar; //cout<<"Pi "<<Pi.rows()<<Pi.cols()<<endl;
	//MatrixXd Eye; Eye.setIdentity(Pi.rows(),Pi.cols());
	MatrixType Eye(MatrixType::Identity(Pi.rows(),Pi.cols())); //cout<<"Eye "<<Eye.rows()<<Eye.cols()<<endl;
	//cout<<"Gtilda "<<ComputeGTilda(k).rows()<<ComputeGTilda(k).cols()<<endl;
	MatrixType ret=Pistar.transpose()*ComputeGTilda(k)*Pistar+(Eye-Pi).transpose()*(Eye-Pi);
	cout<<"Local stiffness dimensions"<<ret.rows()<<" "<<ret.cols()<<endl;
	return ret;
}