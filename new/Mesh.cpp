#include "Mesh.hpp"
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

/*
MeshHandler::MeshHandler(Mesh &mesh):
pointList(mesh.M_pointList),
elementList(mesh.M_elementList),
edgeList(mesh.M_edgeList),
bEdgeList(mesh.M_bEdgeList),
m(mesh)
{};
*/


/*
int MeshRead::read(Mesh & m, std::string const & filename){
	using namespace std;
	MeshHandler mesh(m);
	vector<Point> & pl(mesh.pointList);
	vector<Triangle> & el(mesh.elementList);
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
		streampos oldpos=file.tellg();
		getline(file,currLine);
		stringstream ss(currLine);
		double X,Y; ss>>X>>Y;
		char cha; ss>>cha;
		//ora se tutto va bene sono alla fine della stringa. Se non lo sono ho finito
		if (!ss.eof()) {cout<<"Arrived at number "<<M_pointList.size()<<endl; file.seekg(oldpos); break;}

		//aggiungi il punto trovato
		M_pointList.push_back(Point2D{X,Y});

		if(this->M_verbose) 
			cout<<"Added the point number "<<M_pointList.size()<<" which is X= "<<X<<" Y= "<<Y<<endl;
	}

	//reads connectivity matrix
	while(1){
		streampos oldpos=file.tellg();
		getline(file,currLine);
		stringstream ss(currLine);

		//salva la i-esima riga della matrice
		std::vector<unsigned int> line;
		unsigned int d;
		while (!ss.eof()) {ss>>d; line.push_back(d);}
		line.pop_back(); //se non lo metto, l'ultimo elemento viene contato due volte

		//se tutto va bene, la lunghezza deve essere almeno 2. Se non lo è ho finito
		if (line.size()<=1) {cout<<"Arrived at number "<<M_elementList.size()<<endl; file.seekg(oldpos); break;}


		for (unsigned int i=1; i<line.size(); ++i){
			M_edgeList.insert(Edge{line[i],line[i-1]});
			if(this->M_verbose) 
				cout<<"Added the edge "<<M_edgeList.size()<<" which is a= "<<line[i]<<" b= "<<line[i-1]<<endl;
		}
		M_edgeList.insert(Edge{line[0],line[line.size()]});
		if(this->M_verbose) 
			cout<<"Added the edge "<<M_edgeList.size()<<" which is a= "<<line[0]<<" b= "<<line[line.size()]<<endl;

		M_elementList.push_back(Polygon(line));
		if(this->M_verbose) 
			cout<<"Added the polygon "<<M_elementList.size()<<" which is "<<Polygon(line).showMe()<<endl;
	}

	//reads boundary elements (dovrei salvare gli edges, ma qui ho solo le coordinate dei vertici)
	//boundary dovrebbe essere un vettore di unsigned int
	while (getline(file,currLine)){
		boundary.push_back(stoi(currLine));
		if(this->M_verbose) 
			cout<<"New boundary vertex "<<boundary.size()<<" which is "<<stoi(currLine)<<endl;
	}


};


  
Mesh::Mesh(std::string const filename, MeshReader & reader)
{reader.read(*this,filename);}

  

int Mesh::readMesh(const std::string & file, MeshReader & reader)
{return reader.read(*this,file);}

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
  double Mesh::measure()const 
  {
    double mes(0);
    for (auto i=M_elementList.cbegin();
	 i!=M_elementList.cend();++i)mes+=i->measure();
    return mes;
  }

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




