#include "grid.hpp"
#include <string>
#include <algorithm>

using namespace Geometry;
using namespace std;

//    std::vector<Point2D> coord;
//    std::vector<std::shared_ptr<AbstractPolygon> > abspol;
//    std::vector<unsigned int> boundary;

Grid::Grid(ifstream & file){
string str;

//reads coordinates of the vertexes
while (1){
	streampos oldpos=file.tellg();
	getline(file,str);
	stringstream ss(str);
	double X,Y;
	ss>>X>>Y;
	char cha;
	ss>>cha;
	if (!ss.eof()) {cout<<"Arrived at number "<<coord.size()<<endl; file.seekg(oldpos); break;}
	coord.push_back(Point2D{X,Y});
	//cout<<coord.size()<<endl;
}

//reads connectivity matrix
while(1){
	streampos oldpos=file.tellg();
	getline(file,str);
	stringstream ss(str);
	std::vector<unsigned int> line;
	unsigned int d;
	//char cha;
	while (!ss.eof()) {ss>>d; line.push_back(d);}
	line.pop_back();

	if (line.size()<=1) {cout<<"Arrived at number "<<abspol.size()<<endl; file.seekg(oldpos); break;}

	for (unsigned int i=1; i<line.size(); i++){
		//cout<<"Inserisco "<<line[i]<<""<<line[i-1]<<endl;
		edges.insert(Edge{line[i],line[i-1]});
	}
	edges.insert(Edge{line[0],line[line.size()]});

	vector<Point2D> p;
	for (unsigned int j=0; j<line.size(); j++) {p.push_back(Point2D(coord[line[j]-1]));}
	abspol.push_back(make_shared<Polygon>(Polygon(p)));
	//Polygon(p).showMe();
	//cout<<"Arrived at polygon"<<abspol.size()<<endl;
	//cout<<line[0]<<" "<<line[1]<<" "<<line[2]<<line[3]<<endl;
}

//reads boundary elements
//int i=0;
while (getline(file,str)){
	boundary.push_back(stoi(str));
	//cout<<boundary[i]<<endl; i++;
}

for (auto i : edges) {cout<<i;}
}







double Grid::area(){
	unsigned int Npoly=abspol.size();
	double area{0.0};
	for (unsigned int i=0; i<Npoly; i++){
		area+=abspol[i]->area();
	}
return area;
}

void Grid::ConnMatrixBound(int k){
	for (unsigned int i=0; i<abspol.size(); i++){
		Vertices v=abspol[i]->BoundaryDof(k);
		//for (auto ii : v) cout<<ii.x()<<" "<<ii.y()<<endl;
		for (auto j : v) connBound.insert(j);
		//cout<<"End Polygon"<<endl;
		//cout<<"Actual size"<<connBound.size()<<endl;
	}
}

MatrixType Grid::K(int k){
	//dimensione: interni+num vertici+ boundary
	unsigned int dim=abspol.size()*(k-1)*(k)/2+coord.size()+(k-1)*edges.size();
	MatrixType K(dim,dim);
	cout<<"Dimensions "<<dim<<endl;
	ConnMatrixBound(k);
	//for (auto i : connBound) cout<<i.x()<<" "<<i.y()<<endl;
	cout<<"boundary dofs: "<<connBound.size()<<endl;


	//for (unsigned int cont=0; cont<abspol.size(); cont++){
	for (unsigned int cont=0; cont<abspol.size(); cont++){
		MatrixType locK=(abspol[cont])->ComputeStiffness(k);
			cout<<" Arrived here"<<endl;

		//assemble line of connettivity matrix (parte già da 0!!!)
		vector<unsigned int> line(((*(abspol[cont])).theVertices()).size());
		//cout<<"Line size"<<line.size()<<endl;
		
			for (unsigned int i=0; i<line.size(); i++){
				Point2D PP(((abspol[cont])->theVertices())[i]);
				auto it=find(coord.begin(),coord.end(),PP);
				line[i]=distance(coord.begin(),it);
			}

		cout<<"Current line: "<<line[0]<<" "<<line[1]<<" "<<line[2]<<" "<<line[3]<<endl;

		//assemble line of connectivity matrix for boundary dof (check)
		vector<unsigned int> line2;
		Vertices v=abspol[cont]->BoundaryDof(k);
		for (unsigned int j=0; j<v.size(); j++) {
			auto it=find(connBound.begin(),connBound.end(),v[j]);
			line2.push_back(distance(connBound.begin(),it));
		}
		cout<<"Current line: "<<line2[0]<<" "<<line2[1]<<" "<<line2[2]<<" "<<line2[3]<<endl;

		//internal dofs (to be checked)
		vector<unsigned int> line3;
		for (unsigned int i=0; i<k*(k-1)/2; i++) line3.push_back(cont+i);

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

