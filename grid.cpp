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

	for (unsigned int i=1; i<line.size(); i++)
		edges.insert(Edge{line[i],line[i-1]});

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

}







double Grid::area(){
	unsigned int Npoly=abspol.size();
	double area{0.0};
	for (unsigned int i=0; i<Npoly; i++){
		area+=abspol[i]->area();
	}
return area;
}



MatrixType Grid::K(int k){
	unsigned int dim=abspol.size()*k*(k-1)/2+coord.size()+(k-1)*edges.size();
	MatrixType K(dim,dim);

/*
	for (unsigned int cont=0; cont<abspol.size(); cont++){
		MatrixType locK=(abspol[cont])->ComputeStiffness(k);
			cout<<" Arrived here"<<endl;

		//assemble line of connettivity matrix
		vector<unsigned int> line(((*(abspol[cont])).theVertices()).size());
		//cout<<"Line size"<<line.size()<<endl;
		
			for (unsigned int i=0; i<line.size(); i++){
				Point2D PP(((abspol[cont])->theVertices())[i]);
				auto it=find(coord.begin(),coord.end(),PP);
				line[i]=distance(coord.begin(),it);
			}

		cout<<"Current line: "<<line[0]<<" "<<line[1]<<" "<<line[2]<<line[3]<<endl;
		

		for (unsigned int i=0; i<locK.rows(); i++){
			for (unsigned int j=0; j<locK.cols(); j++){
				//unsigned int index1=((*(abspol[cont])).theVertices())[i]
				K(line[i],line[j])
			}
		}

	}
*/
	return K;
};

