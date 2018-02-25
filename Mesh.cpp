#include "Mesh.hpp"
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <set>
#include "freefunc.hpp"


MeshHandler::MeshHandler(Mesh &mesh):
pointList(mesh.M_pointList),
elementList(mesh.M_elementList),
boundary(mesh.M_boundary),
num_edges(mesh.M_num_edges),
m(mesh)
{};




int MeshReader::read(Mesh & m, std::string const & filename){
	using namespace std;
	std::vector<unsigned int> tmp;
	MeshHandler mesh(m);
	vector<Point> & pl(mesh.pointList);
	vector<Polygon> & el(mesh.elementList);
	ifstream f;
	string currLine;

	f.open(filename.c_str());
	if (!f.is_open()) {
		cerr<<"Mesh file does not exist or is corrupted"<<endl;
		return 1;
	}

	if(this->M_verbose) cout<<"Reading mesh data"<<endl;

	//reads coordinates of the vertexes
	while (1){
		streampos oldpos=f.tellg();
		getline(f,currLine);
		stringstream ss(currLine);
		double X,Y; ss>>X>>Y;
		char cha; ss>>cha;

		//if I am ok I reached end of the string. otherwise go back to beginning of line and break
		if (!ss.eof()) {cout<<"Total number of points = "<<pl.size()<<endl; f.seekg(oldpos); break;}

		//add point
		pl.push_back(Point{X,Y});

		if(this->M_verbose) 
			cout<<"Added the point number "<<pl.size()-1<<" which is X= "<<X<<" Y= "<<Y<<endl;
	}

	//reads connectivity matrix
	while(1){
		streampos oldpos=f.tellg();
		getline(f,currLine);
		stringstream ss(currLine);

		//save i-th line of the matrix
		//note: connectivity matrix in input starts by 1, but I need it to start from 0
		std::vector<unsigned int> line;
		unsigned int d;
		while (!ss.eof()) {ss>>d; line.push_back(d-1);}
		line.pop_back(); //otherwise I count the last element twice

		//if I am ok, the length has to be at least 2. otherwise go back to beginning of line and break
		if (line.size()<=1) {cout<<"Total number of polygons = "<<el.size()<<endl; f.seekg(oldpos); break;}

		el.push_back(Polygon(line,&pl));

		if(this->M_verbose) 
			cout<<"Added the polygon "<<el.size()-1<<" which is "<<Polygon(line,&pl)<<endl;

	}

	//reads boundary elements (indexes, with shift by one)
	while (getline(f,currLine)){
		
		mesh.boundary.push_back(stoi(currLine)-1);
		
		if(this->M_verbose) 
			cout<<"Added boundary vertex with index "<<mesh.boundary.size()-1<<" which is "<<stoi(currLine)-1<<endl;
	}
	cout<<"Total number of boundary vertexes = "<<mesh.boundary.size()<<endl;

	return 0;

};


  
Mesh::Mesh(std::string const filename, MeshReader & reader, unsigned int kk)
{k=kk; reader.read(*this,filename);}

  

int Mesh::readMesh(const std::string & file, MeshReader & reader, unsigned int kk)
{k=kk; return reader.read(*this,file);}


std::ostream & operator << (std::ostream & ost, const Mesh & m){
	ost<<"##MESH##"<<std::endl;
	ost<<"#POINTS#"<<std::endl;
	for (auto i : m.M_pointList) ost<<i;
	ost<<"#POLYGONS#"<<std::endl;
	for (auto i : m.M_elementList) ost<<i;
	ost<<"#BOUNDARY VERTEXES#"<<std::endl;
	for (unsigned int i=0; i<m.M_boundary.size(); ++i) 
		ost<<m.M_boundary[i]<<" corresponding to point "<<m.M_pointList[m.M_boundary[i]];
return ost;
}

double Mesh::area()const {
	double meas(0);
	for (auto i : M_elementList) meas+=i.area();
	return meas;
}

double Mesh::max_diam()const {
	double d(0);
	for (auto i : M_elementList) d=std::max(i.diameter(),d);
	return d;
}




//boundary dofs sets the indexes for the dofs on the edges, not necessarily Dirichlet ones
void Mesh::boundaryDOF(){

	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		auto points=M_elementList[i].getPoints();
		std::vector<Point> nodes;
		std::vector<double> weights;
		computeDOF(points,k,weights,nodes);
		
		std::vector<unsigned int> line;
		for (unsigned int j=0; j<nodes.size(); j++) {
			
			//if is not in the vector, insert it, otherwise gives the index
			auto index=std::find(M_edgesDOF.begin(),M_edgesDOF.end(),nodes[j]);
			if (index==M_edgesDOF.end()) {
				M_edgesDOF.push_back(nodes[j]); line.push_back(M_edgesDOF.size()-1); 
			}
			else 
				line.push_back(std::distance(M_edgesDOF.begin(),index));
				
			//setting dofs for the polygon
			M_elementList[i].setDof(line,&M_edgesDOF);
		}
	}
	return;
}


//diffusion+transport term
MatrixType_S Mesh::GlobalStiffness(std::function<double (double, double)> mu, double mu_bar, bool constant_mu, 
	std::function<double (double, double)> beta_x, std::function<double (double, double)> beta_y){
	
	boundaryDOF(); 
	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType_S K(dim,dim); 
	//K.fill(0.0);
	std::cout<<"Dimension of the global stiffness matrix is "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		MatrixType locK=M_elementList[i].LocalStiffness_weighted(k,mu,mu_bar,constant_mu)+
			M_elementList[i].LocalTransport(k, beta_x, beta_y);
		
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		//std::cout<<"Current dofs: ";
		//for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global stiffness matrix
		for (unsigned int a=0; a<locK.rows(); a++){
			for (unsigned int b=0; b<locK.cols(); b++){
				K.coeffRef(current[a],current[b])+=locK(a,b);
			}
		}
		
	} //end loop i
	
	return K;
}



MatrixType_S Mesh::GlobalMass(){
	boundaryDOF();

	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType_S M(dim,dim); 

	std::cout<<"Dimension of the global mass matrix is "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		MatrixType locM=M_elementList[i].LocalMass(k);
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		//std::cout<<"Current dofs: ";
		//for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global matrix
		for (unsigned int a=0; a<locM.rows(); a++){
			for (unsigned int b=0; b<locM.cols(); b++){
				M.coeffRef(current[a],current[b])+=locM(a,b);
			}
		}
		
		} //end loop i
		
	return M;
}







MatrixType Mesh::GlobalLoad(std::function<double(double,double)> f){
	boundaryDOF();

	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType F(dim,1); 
	F.fill(0.0);
	std::cout<<"Dimension of the global load term is "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		MatrixType locF=M_elementList[i].LoadTerm(k,f); //compute local load term
		
		//find the indexes of the dofs of the polygon w.r.t. global problem
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
	
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		//std::cout<<"Current dofs: ";
		//for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global vector
		for (unsigned int a=0; a<locF.rows(); a++){
			F(current[a],0)+=locF(a,0);
		}
		
		} //end loop i
		
	return F;
}




MatrixType Mesh::solve(std::function<double (double,double)> f, std::function<double (double,double)> g,
	std::function<double (double,double)> mu, double mu_bar, bool constant_mu, 
	std::function<double (double,double)> beta_x, std::function<double (double,double)> beta_y){
		
	//compute global LHS and RHS
	std::cout<<"Computing global stiffness matrix"<<std::endl;
	MatrixType_S K=GlobalStiffness(mu,mu_bar,constant_mu, beta_x, beta_y);
	std::cout<<"Computing global load term"<<std::endl;
	MatrixType F=GlobalLoad(f);
	//compute indexes associated with dofs on the domain boundary to apply Dirichlet BC
	std::vector<unsigned int> Dir=Dirichlet();
	

	int ii=-1, jj=0,jjj=0;

	//final solution
	MatrixType U(K.rows(),1);
	
	//create a boundary solution to store Dirichlet boundary condition
	MatrixType UB(Dir.size(),1);
	ii=0;
	for (unsigned int i=0; i<U.rows(); i++){
		if (find(Dir.begin(),Dir.end(),i)!=Dir.end()) {
			
			//if Dirichlet node, get the point
			Point PP;
			if(i<M_pointList.size()) 
				PP=M_pointList[i];
			else 
				PP=M_edgesDOF[i-M_pointList.size()];
				
			//compute boundary data in that point
			UB(ii,0)=g(PP[0],PP[1]); 
			ii++;
		
		}
	} //end loop i


	//set to zeros entries of the matrix whenever the row index is boundary
	for (int k=0; k<K.outerSize(); ++k)
		for (MatrixType_S::InnerIterator it(K,k); it; ++it) {
		if (find(Dir.begin(),Dir.end(),it.row())!=Dir.end()){
			K.coeffRef(it.row(),it.col())=0.0;
			}
		}

	ii=0;
	for (unsigned int i=0; i<F.rows(); i++){
			if (find(Dir.begin(),Dir.end(),i)!=Dir.end()) {
				K.coeffRef(i,i)=1.0; //identity matrix at boundary indexes
				F(i,0)=UB(ii,0); //the source term is equal to the boundary value
				ii++; 
			}
	}

	//solve the system
	std::cout<<"Solving the linear system"<<std::endl;
	
	//create a vector to store the result of system solution
	MatrixType UI(U.rows()-UB.rows(),1);
	//solve
	Eigen::BiCGSTAB<MatrixType_S> solver; solver.compute(K); U=solver.solve(F);
	
	//gather contributions from UI and UB
	ii=0; unsigned int iii=0;
	for (unsigned int i=0; i<U.rows(); i++){
		
		if (find(Dir.begin(),Dir.end(),i)!=Dir.end()) {
			U(i,0)=UB(ii,0); //if it is Dirichlet take from UB 
			ii++;
		}
	
	} //end loop i
	
	//print VEM solution on screen
	//std::cout<<"Solution:"<<std::endl<<U<<std::endl;

	//print VEM solution on file (usually I need only boundary dofs, integrals are not needed)
	std::ofstream file("output.dat");
	for (unsigned int i=0; i<M_pointList.size(); i++) 
		file<<M_pointList[i][0]<<"\t"<<M_pointList[i][1]<<"\t"<<U(i,0)<<std::endl;
	for (unsigned int i=0; i<M_edgesDOF.size(); i++) 
		file<<M_edgesDOF[i][0]<<"\t"<<M_edgesDOF[i][1]<<"\t"<<U(i+M_pointList.size(),0)<<std::endl;
	//internal dofs (if needed)
	//for (unsigned int i=0; i<k*(k-1)/2*M_elementList.size(); i++) 
	//	file<<"Polygon "<<i/(k*(k-1)/2)<<" dof number "<<i%(k*(k-1)/2)<<" "<<U(i+M_pointList.size()+M_edgesDOF.size(),0)<<std::endl;
	file<<std::endl;
	
	file.close();


	return U;
}




std::vector<unsigned int> Mesh::Dirichlet(){
	
	std::vector<unsigned int> Dir=M_boundary;
	
	//for each polygon, loop over dofs and find if it is a boundary dof
	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		
		//strategy: if both vertexes are boundary dofs, then the all the intermediate GL points are boundary
		//note: this works only if I have at least two elements in each direction (it has to be modified)
		
		for (unsigned int j=0; j<line1.size(); j++){
			//se entrambi gli estremi sono di bordo, allora il dof è di bordo (NON è VERO!)
			if (find(M_boundary.begin(), M_boundary.end(),line1[j])!=M_boundary.end() &&
				find(M_boundary.begin(), M_boundary.end(),line1[(j+1)%line1.size()])!=M_boundary.end())
		
				for (unsigned int z=0; z<k-1; z++) 
					Dir.push_back(line2[j*(k-1)+z]+M_pointList.size());
					 
		}
	
	} //end loop i
	
	//for (auto i : Dir) std::cout<<"Boundary indexes: "<<i<<std::endl;
	
	return Dir;
}


MatrixType Mesh::VEMConvert(std::function<double (double,double)> uex){

	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType U(dim,1); 

	for (unsigned int i=0; i<M_elementList.size(); i++){
		
		MatrixType locU=M_elementList[i].LocalConvert(k,uex);

		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		//std::cout<<"Current dofs: ";
		//for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global vector (note: I do not do += otherwise I sum multiple contributions from the same dof)
		for (unsigned int a=0; a<locU.rows(); a++){
				U(current[a],0)=locU(a,0);
		}
		
	} //end loop i
	
	return U;
}


double Mesh::normInf(MatrixType & uex,MatrixType & u){
	
	//consider only boundary dofs
	unsigned int maxindex=M_pointList.size()+M_edgesDOF.size();
	
	MatrixType diff=uex-u;
	double norm=0.0;
	for (unsigned int i=0; i<maxindex; i++)
		norm=std::max(norm,std::abs(diff(i,0)));

	return norm;
}



void Mesh::Allnorms(MatrixType & u, MatrixType & uex){
	
	
	MatrixType_S K=GlobalStiffness([](double x, double y) {return 1.0;},1.0,true,
		[](double x, double y) {return 0.0;}, [](double x, double y) {return 0.0;});
	MatrixType_S M=GlobalMass();

	//std::cout<<"Exact solution (converted VEM): "<<std::endl<<uex<<std::endl;
	//std::cout<<"VEM approximation of exact solution (numerical): "<<std::endl<<u<<std::endl;

	std::cout<<"Infinity norm: "<<normInf(uex,u)<<std::endl;

	MatrixType H10=((u-uex).transpose())*K*(u-uex);
	std::cout<<"H1 seminorm: "<<std::sqrt(H10(0,0))<<std::endl;

	MatrixType L2=((u-uex).transpose())*M*(u-uex);
	std::cout<<"L2 norm: "<<std::sqrt(L2(0,0))<<std::endl;

	std::cout<<"H1 norm: "<<std::sqrt(H10(0,0)+L2(0,0))<<std::endl;
	return;
}



