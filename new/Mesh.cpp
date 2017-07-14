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
//edgeList(mesh.M_edgeList),
//bEdgeList(mesh.M_bEdgeList),
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

	//unsigned int & num_edges(mesh.num_edges);

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
		//ora se tutto va bene sono alla fine della stringa. Se non lo sono ho finito
		if (!ss.eof()) {cout<<"Total number of points = "<<pl.size()<<endl; f.seekg(oldpos); break;}

		//aggiungi il punto trovato
		pl.push_back(Point{X,Y});

		if(this->M_verbose) 
			cout<<"Added the point number "<<pl.size()-1<<" which is X= "<<X<<" Y= "<<Y<<endl;
	}

	//reads connectivity matrix
	while(1){
		streampos oldpos=f.tellg();
		getline(f,currLine);
		stringstream ss(currLine);

		//salva la i-esima riga della matrice
		//nota: la matrice di connettività parte da 1, a me serve che parta da 0!!!
		std::vector<unsigned int> line;
		unsigned int d;
		while (!ss.eof()) {ss>>d; line.push_back(d-1);}
		line.pop_back(); //se non lo metto, l'ultimo elemento viene contato due volte

		//se tutto va bene, la lunghezza deve essere almeno 2. Se non lo è ho finito
		if (line.size()<=1) {cout<<"Total number of polygons = "<<el.size()<<endl; f.seekg(oldpos); break;}

		/*
		for (unsigned int i=1; i<line.size(); ++i){
			M_edgeList.insert(Edge{line[i],line[i-1]});
			if(this->M_verbose) 
				cout<<"Added the edge "<<M_edgeList.size()<<" which is a= "<<line[i]<<" b= "<<line[i-1]<<endl;
		}
		M_edgeList.insert(Edge{line[0],line[line.size()]});
		if(this->M_verbose) 
			cout<<"Added the edge "<<M_edgeList.size()<<" which is a= "<<line[0]<<" b= "<<line[line.size()]<<endl;
		*/
		el.push_back(Polygon(line,&pl));
		if(this->M_verbose) 
			cout<<"Added the polygon "<<el.size()-1<<" which is "<<Polygon(line,&pl)<<endl;

		//devo capire se ho già inserito il lato o no (to be done)
		/*
		for (unsigned int j=0; j<line.size(); j++){
			if (find(tmp.begin(),tmp.end(),line[j])==tmp.end()) tmp.push_back(line[j]);
			  find(tmp.begin(),tmp.end(),line[(j+1)%line.size()])==tmp.end()
				{tmp.push_back(line[j]); tmp.push_back(line[(j+1)%line.size()]); num_edges++;}
		}
		*/
	}

	//reads boundary elements (dovrei salvare gli edges, ma qui ho solo le coordinate dei vertici)
	//anche qui shift di 1
	while (getline(f,currLine)){
		mesh.boundary.push_back(stoi(currLine)-1);
		if(this->M_verbose) 
			cout<<"Added boundary vertex with index "<<mesh.boundary.size()-1<<" which is "<<stoi(currLine)-1<<endl;
	}
	cout<<"Total number of boundary vertexes = "<<mesh.boundary.size()<<endl;
	//cout<<"Total number of edges = "<<num_edges<<endl;

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





void Mesh::boundaryDOF(){
	//std::set<Point> M_set;
	//M_edgesDOF.clear();
	for (unsigned int i=0; i<M_elementList.size(); i++){
		auto points=M_elementList[i].getPoints();
		std::vector<Point> nodes;
		std::vector<double> weights;
		computeDOF(points,k,weights,nodes);
		std::vector<unsigned int> line;
		for (unsigned int j=0; j<nodes.size(); j++) {
			//insert if is not in the vector
			auto index=std::find(M_edgesDOF.begin(),M_edgesDOF.end(),nodes[j]);
			if (index==M_edgesDOF.end())
				{M_edgesDOF.push_back(nodes[j]); line.push_back(M_edgesDOF.size()-1); 
					//std::cout<<std::distance(M_edgesDOF.begin(),index+1)<<"hh"<<std::endl;
				}
			else {line.push_back(std::distance(M_edgesDOF.begin(),index));
					//std::cout<<std::distance(M_edgesDOF.begin(),index)<<"hh"<<std::endl;
				}
			//auto status=M_set.insert(nodes[j]);
			//for (auto kk : M_edgesDOF) std::cout<<kk<<std::endl;
			//std::cout<<"Inserted element number "<<line[line.size()-1]<<"corresponding to point "<<nodes[j];
			M_elementList[i].setDof(line,&M_edgesDOF);
		}
	}
	return;
}



MatrixType Mesh::GlobalStiffness(){
	boundaryDOF(); //std::cout<<"Number of BD = "<<M_edgesDOF.size()<<std::endl;
	//dimensione: dof interni + numero vertici + boundary + dof interni
	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType K(dim,dim); K.fill(0.0);
	std::cout<<"Dimensions "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		MatrixType locK=M_elementList[i].LocalStiffness(k);
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		std::cout<<"Current dofs: ";
		for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global matrix
		for (unsigned int a=0; a<locK.rows(); a++){
			for (unsigned int b=0; b<locK.cols(); b++){
				//std::cout<<"Inserisco elemento "<<a<<" "<<b<<" in posizione "<<current[a]<<" "<<current[b]<<std::endl;
				K(current[a],current[b])+=locK(a,b);
			}
		}
		} //end for of polygons
	return K;
}



MatrixType Mesh::GlobalMass(){
	boundaryDOF(); //std::cout<<"Number of BD = "<<M_edgesDOF.size()<<std::endl;

	//dimensione: dof interni + numero vertici + boundary + dof interni
	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType M(dim,dim); M.fill(0.0);
	std::cout<<"Dimensions "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		MatrixType locM=M_elementList[i].LocalMass(k);
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		std::cout<<"Current dofs: ";
		for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global matrix
		for (unsigned int a=0; a<locM.rows(); a++){
			for (unsigned int b=0; b<locM.cols(); b++){
				M(current[a],current[b])+=locM(a,b);
			}
		}
		} //end for of polygons
	return M;
}



MatrixType Mesh::GlobalLoad(std::function<double(double,double)> f){
	boundaryDOF(); //std::cout<<"Number of BD = "<<M_edgesDOF.size()<<std::endl;

	//dimensione: dof interni + numero vertici + boundary + dof interni
	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType F(dim,1); F.fill(0.0);
	std::cout<<"Dimensions "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		MatrixType locF=M_elementList[i].LoadTerm(k,f);
		//std::cout<<locF;
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		std::cout<<"Current dofs: ";
		for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global vector
		for (unsigned int a=0; a<locF.rows(); a++){
			F(current[a],0)+=locF(a,0);
		}
		} //end for of polygons
	return F;
}




MatrixType Mesh::solve(std::function<double (double,double)> f, std::function<double (double,double)> g){
	MatrixType K=GlobalStiffness();
	MatrixType F=GlobalLoad(f);
	std::vector<unsigned int> Dir=Dirichlet();
	MatrixType KII(K.rows()-Dir.size(),K.cols()-Dir.size());
	std::cout<<KII.rows()<<"  "<<KII.cols()<<std::endl;
	MatrixType KIB(K.rows()-Dir.size(),Dir.size());
	std::cout<<KIB.rows()<<"  "<<KIB.cols()<<std::endl;
	MatrixType FI(F.rows()-Dir.size(),1);
	int ii=-1, jj=0,jjj=0;

	for (unsigned int i=0; i<K.rows(); i++){
		//std::cout<<i<<std::endl;
		jj=0; jjj=0;
		if (find(Dir.begin(),Dir.end(),i)==Dir.end()) {
		ii++; //se la i non è di Dirichlet aumenta di 1
		for (unsigned int j=0; j<K.cols(); j++){
			if (find(Dir.begin(),Dir.end(),j)==Dir.end()) //se anche j non è di Dirichlet inserisci in KII
				{KII(ii,jj)=K(i,j); //std::cout<<"internal"<<ii<<"  "<<jj<<std::endl; 
				jj++;
				}
			else {KIB(ii,jjj)=K(i,j); //std::cout<<"boundary"<<ii<<"  "<<jjj<<std::endl; 
				jjj++;
					}
		}
		}
	} //fine for

	//termine noto
	ii=0;

	std::cout<<"Rows "<<F.rows();
	for (unsigned int i=0; i<F.rows(); i++){
		//std::cout<<"hey"<<i<<std::endl;
		if (find(Dir.begin(),Dir.end(),i)==Dir.end()) {FI(ii,0)=F(i,0); ii++;}
	}
	//std::cout<<FI<<std::endl;

	//real solution
	MatrixType U(K.rows(),1), UB(Dir.size(),1);
	ii=0;
	for (unsigned int i=0; i<U.rows(); i++){
		//std::cout<<i<<std::endl;
		if (find(Dir.begin(),Dir.end(),i)!=Dir.end()) {
			Point PP;
			if(i<M_pointList.size()) PP=M_pointList[i];
			else PP=M_edgesDOF[i-M_pointList.size()];
			UB(ii,0)=g(PP[0],PP[1]); 
			ii++;}
	}
	//std::cout<<UB<<std::endl;

	//solve
	MatrixType UI(U.rows()-UB.rows(),1);
	UI=(KII.lu()).solve(FI-KIB*UB);
	//assemble back
	ii=0; unsigned int iii=0;
	for (unsigned int i=0; i<U.rows(); i++){
		//std::cout<<i<<std::endl;
		if (find(Dir.begin(),Dir.end(),i)!=Dir.end()) {
			U(i,0)=UB(ii,0); ii++;}
		else {U(i,0)=UI(iii,0); iii++;}
	}
	std::cout<<"Solution:"<<std::endl<<U<<std::endl;
	//std::cout<<"MATRICE"<<KIB;

	//print solution
	std::ofstream file("output.dat");
	for (unsigned int i=0; i<M_pointList.size(); i++) 
		file<<M_pointList[i][0]<<"\t"<<M_pointList[i][1]<<"\t"<<U(i,0)<<std::endl;
	for (unsigned int i=0; i<M_edgesDOF.size(); i++) 
		file<<M_edgesDOF[i][0]<<"\t"<<M_edgesDOF[i][1]<<"\t"<<U(i+M_pointList.size(),0)<<std::endl;
	//internal DOFs
	//for (unsigned int i=0; i<k*(k-1)/2*M_elementList.size(); i++) 
	//	file<<"Polygon "<<i/(k*(k-1)/2)<<" dof number "<<i%(k*(k-1)/2)<<" "<<U(i+M_pointList.size()+M_edgesDOF.size(),0)<<std::endl;
	file<<std::endl;
//end function
return U;
}

std::vector<unsigned int> Mesh::Dirichlet(){
	std::vector<unsigned int> Dir=M_boundary;
	for (unsigned int i=0; i<M_elementList.size(); i++){
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		for (unsigned int j=0; j<line1.size(); j++){
			//se entrambi gli estremi sono di bordo, allora il dof è di bordo (NON è VERO!)
			if (find(M_boundary.begin(), M_boundary.end(),line1[j])!=M_boundary.end() &&
				find(M_boundary.begin(), M_boundary.end(),line1[(j+1)%line1.size()])!=M_boundary.end())
				{for (unsigned int z=0; z<k-1; z++) Dir.push_back(line2[j*(k-1)+z]+M_pointList.size());} 
		}
	}
	for (auto i : Dir) std::cout<<"Indici di bordo "<<i<<std::endl;
	return Dir;
}


MatrixType Mesh::VEMConvert(std::function<double (double,double)> uex){

	//costruisco vettore che approssima con i vem
	unsigned int dim=M_pointList.size()+M_edgesDOF.size()+M_elementList.size()*(k-1)*(k)/2;
	MatrixType U(dim,1); U.fill(0.0);
	std::cout<<"Dimensions "<<dim<<std::endl;

	for (unsigned int i=0; i<M_elementList.size(); i++){
		MatrixType locU=M_elementList[i].LocalConvert(k,uex);
		std::cout<<locU;
		std::vector<unsigned int> line1=M_elementList[i].getVertexes();
		std::vector<unsigned int> line2=M_elementList[i].getBDindexes();
		std::vector<unsigned int> line3;
		for (unsigned int j=0; j<k*(k-1)/2; j++) line3.push_back(j); 
		std::vector<unsigned int> current=line1;
		for (unsigned int j=0; j<line2.size(); j++) current.push_back(line2[j]+M_pointList.size());
		for (unsigned int j=0; j<line3.size(); j++) current.push_back(line3[j]+M_pointList.size()+M_edgesDOF.size()+i*k*(k-1)/2);

		std::cout<<"Current dofs: ";
		for (auto j : current) std::cout<<j<<" "; std::cout<<std::endl;

		//assemble global vector (mi basta un singolo integrale, l'if è indifferente metterlo o no)
		for (unsigned int a=0; a<locU.rows(); a++){
			//if (current[a]>=line1.size()+line2.size()) ????
			//	U(current[a],0)+=locU(a,0);
			//else 
				U(current[a],0)=locU(a,0);
		}
		} //end for of polygons
	return U;
}


double Mesh::normInf(MatrixType uex,MatrixType u){
	unsigned int maxindex=M_pointList.size()+M_edgesDOF.size();
	MatrixType diff=uex-u;
	double norm=0.0;
	for (unsigned int i=0; i<maxindex; i++)
		norm=std::max(norm,std::abs(diff(i,0)));
	//std::cout<<"Expected "<<norm<<std::endl;
	/*
	maxindex=diff.rows(); norm=0.0;
	for (unsigned int i=0; i<maxindex; i++)
		norm=std::max(norm,std::abs(diff(i,0)));
	std::cout<<"wrong"<<norm<<std::endl;
	*/
	return norm;
}

double Mesh::H1seminorm(MatrixType uex, MatrixType u, MatrixType K){
	MatrixType temp=((u-uex).transpose())*K*(u-uex);
	return std::sqrt(temp(0,0));
}


void Mesh::Allnorms(MatrixType u, MatrixType uex){
	MatrixType K=GlobalStiffness(), M=GlobalMass();

	std::cout<<"Exact solution (converted VEM): "<<std::endl<<uex<<std::endl;
	std::cout<<"VEM approximation of exact solution (numerical): "<<std::endl<<u<<std::endl;

	std::cout<<"Infinity norm: "<<normInf(uex,u)<<std::endl;

	MatrixType tmp=((u-uex).transpose())*K*(u-uex);
	std::cout<<"H1 seminorm: "<<std::sqrt(tmp(0,0))<<std::endl;

	MatrixType temp=((u-uex).transpose())*M*(u-uex);
	std::cout<<"L2 norm: "<<std::sqrt(temp(0,0))<<std::endl;
	//std::cout<<tmp(0,0)<<" "<<temp(0,0)<<std::endl;
	std::cout<<"H1 norm: "<<std::sqrt(tmp(0,0)+temp(0,0))<<std::endl;
	return;
}

/*


		vector<unsigned int> V=line;
		for (auto i : line2) V.push_back(i+coord.size());
		for (auto i : line3) V.push_back(i+coord.size()+(k-1)*edges.size());
		for (auto i : V) cout<<i<<" ";
			cout<<endl;
		
		for (unsigned int i=0; i<locK.rows(); i++){
			for (unsigned int j=0; j<locK.cols(); j++){
				K(V[i],V[j])+=locK(i,j);
			}
		}
		

	}

	return K;
};
*/


/*
	// Scan file lines
    // Skip until data or end of file
      skipInput(f,std::string("DATA"));
      // If the string is not found, abort
      if(f.eof()||f.fail()){
	std::cerr << "FILE ERROR! Cannot find #DATA" << std::endl;
	return 2;
      }
      // Read number of points and elements
      int nP, nE;
      f >> nP >> nE;
      if(this->M_verbose)std::cout<<"Number of points="<<nP<<
			   " Number of elements="<<nE<<std::endl;
      pl.reserve(nP);
      el.reserve(nE);
      // Now look for the points
      skipInput(f,std::string("POINTS"));
      // If the string is not found, abort
      if(f.eof()||f.fail()||f.bad())
      	{ std::cerr << "FILE ERROR! cannot find # POINTS" << std::endl;
	  return 3; }
      
      // Fill point data structures
      for(int i = 0; i < nP; ++i) {
      	Point P;
      	double x,y;
      	int id;
      	f >> x>>y>>id;
      	P[0]=x;
      	P[1]=y;
      	P.bcId()=id;
      	//std::cout<<x<<" " <<y<<" "<< id<<" "<<pl.size()<<std::endl;
      	//std::cout<<P[0]<<" " <<P[1]<<" "<< P.bcId()<<" "<<pl.size()<<std::endl;
      	P.id()=pl.size();
      	pl.push_back(P);
	if(f.eof()||f.fail())
	  { std::cerr << "FILE ERROR! cannot read all point" << std::endl;
	    return 3; }
      }
      if(this->M_verbose){
	std::cout<<"Points read"<<std::endl;
	std::cout << "Reading elements" << std::endl;
      }
      skipInput(f,std::string("ELEMENTS"));
      
      // If the string is not found, abort
      if(f.eof()||f.fail())
      	{ std::cerr << "FILE ERROR! Cannot find # ELEMENTS" << std::endl;
	  return 4; }
      // Fill element data structures
      for(int i = 0; i < nE; i++) {
      	int type;
        f >> type;
        if(type != 0) {
	  std::cerr << "I can read only triangular elements"<< std::endl;
	  return 5;
        }
        int i1,i2,i3;
        if (!f.eof())f >> i1 >> i2 >> i3;
        if(f.eof()||f.fail()){
	  std::cerr << "FILE ERROR! Cannot read elements" << std::endl;
	  return 5; }
        Triangle t(pl[i1],pl[i2],pl[i3]);
        t.id()=el.size();
        t.bcId()=0;
        el.emplace_back(std::move(t));
      }
      f.close();
      return 0;
      */



/*


  bool Mesh::checkmesh() const{
    using std::clog;
    using std::endl;
    using std::cerr;
    bool status(false);
    clog<<"Mesh stores: "<<endl;
    clog<<this->num_elements()<<" Elements"<<endl;
    clog<<this->num_points()<<" Points"<<endl;
    clog<<this->num_edges()<<" Edges"<<endl;
    clog<<this->num_bEdges()<<" Boundary Edges"<<endl;
    int count(0);
    int bccount(0);
    int wrongid(0);
    size_type idcount(0);
    for(std::vector<Point>::const_iterator i=M_pointList.begin();
	i!=M_pointList.end(); ++i)
      {
	if (i->unassignedId())++count;
	if (i->unassignedBc())++bccount;
	if (i->id()!=idcount++)++wrongid;
	
      }
    clog<<"Mesh has "<<count<<" unassigned point id and "<<
      bccount<<" unassigned point bc markers"<<endl;
    status=wrongid>0;
    if(wrongid>0)cerr<<"Mesh has "<<wrongid<<" wrong point id set";
    count=0;
    bccount=0;
    wrongid=0;
    idcount=0;
    int punset(0);
    for(std::vector<Triangle>::const_iterator i=M_elementList.begin();
	i!=M_elementList.end(); ++i)
      {
	if (i->unassignedId())++count;
	if (i->unassignedBc())++bccount;
	if (i->id()!=idcount++)++wrongid;
	if (i->empty())++punset;Tria
      }
    clog<<"Mesh has "<<count<<" unassigned element id and "<<
      bccount<<" unassigned element bc markers"<<endl;
    status=status||wrongid>0||punset>0;
    if(wrongid>0)cerr<<"Mesh has "<<wrongid<<" wrong element id set";
    if(punset>0)cerr<<"Mesh has "<<punset<<" points unset";
    if(!status)clog<<"Domain area:"<<this->measure()<<endl;
    return status;
  }

  std::ostream & operator<<(std::ostream & out, Mesh const& m)
  {
    out<< " ***** MESH  INFORMATION ******"<<std::endl;
    out<<" Num Points="<<m.num_points()<<" "<<" Num elements="<<m.num_elements()<<" "
       <<"Num. edges="<<m.num_edges()<<" "<<"Num Boundary Edges="
       <<m.num_bEdges()<<std::endl;
    out<< "POINTS:"<<std::endl;
    int oprec=out.precision(10);
    std::ios_base::fmtflags oflags=
      out.setf(std::ios_base::scientific,std::ios_base::floatfield);
    for (std::size_t i=0;i<m.num_points();++i){
      Point p=m.point(i);
      double x=p[0];
      double y=p[1];
      out<<i<<" "<<x<<" "<<y<<std::endl;
    }
    out<<" TRIANGLE CONNECTIVITY AND AREA:"<<std::endl;
    for (std::size_t i=0; i<m.num_elements();++i){
      Triangle t=m.triangle(i);
      out<<i<<" "<<t[0].id()<<" "<<t[1].id()<<" "<<t[2].id()<<
	" "<<t.measure()<<std::endl;
    }
    out.precision(oprec);
    out.flags(oflags);
    return out;
  }
*/




