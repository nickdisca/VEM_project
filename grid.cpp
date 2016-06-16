#include "grid.hpp"
#include <string>

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

	vector<Point2D> p;
	for (unsigned int j=0; j<line.size(); j++) p.push_back(Point2D(coord[line[j]-1]));
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

/*
//coordinates of the vertexes
for (unsigned int i=0; i<Nvert; i++){
	getline(file,str);
	stringstream ss(str);
	ss>>i1>>i1>>i2;
	//cout<<i1<<"space"<<i2<<endl;
	//vect[i].set(stod(aux1),stod(aux2));
	vect[i].set(i1,i2);
}

//stampa
//for (auto i = vect.begin(); i != vect.end(); i++) cout<<i->x()<<" "<<i->y()<<endl;

//connectivity matrix and polygons
for (unsigned int i=0; i<Npoly; i++){
	getline(file,str);
	stringstream ss(str);
	getline(ss,aux,' ');
	getline(ss,aux,' ');
	tipo=stoi(aux);
	vvv.clear(); lati.clear();
	while (getline(ss,aux,' ')){
		vvv.emplace_back(vect[stoi(aux)]);
		lati.emplace_back(stoi(aux));
		//cout<<stoi(aux)<<endl;
	}

	//create a set containing the edges and an auxiliary one to store internal edges
	for(unsigned int j=0; j<lati.size()-1; j++){
		e.set(lati[j],lati[j+1]);
		check=insieme.insert(e);
		if (!check.second) interni.insert(e);
	}
	e.set(lati[lati.size()-1],lati[0]);
	check=insieme.insert(e);
	if (!check.second) interni.insert(e);

	//build the vector of smart pointers
	switch(tipo){
		//case(0): {Triangle t(vvv); abspol.emplace_back(make_shared<Triangle> (t));} break;
		case(0): {abspol.emplace_back(make_shared<Triangle> (Triangle(vvv)));} break;
		case(1): {abspol.emplace_back(make_shared<Square> (Square(vvv)));} break;
		default: {abspol.emplace_back(make_shared<Polygon> (Polygon(vvv)));}
	}
}

//copy all edges in the vector
AllEdges.resize(insieme.size());
copy(insieme.begin(),insieme.end(),AllEdges.begin());

//create the vector containing Boundary edges
Boundary.resize(insieme.size()-interni.size());
set_difference(insieme.begin(), insieme.end(),interni.begin(),interni.end(),Boundary.begin());
//set_difference(insieme.begin(), insieme.end(),interni.begin(),interni.end(),back_inserter(Boundary));
*/
}







double Grid::area(){
	unsigned int Npoly=abspol.size();
	double area{0.0};
	for (unsigned int i=0; i<Npoly; i++){
		area+=abspol[i]->area();
	}
return area;
}


/*
void Grid::printedges(){
	cout<<"All Edges:"<<endl;
	for(auto i=AllEdges.begin(); i!=AllEdges.end(); i++)
		cout<<*i;
	
	cout<<"Lati di bordo:"<<endl;
	for(auto i=Boundary.begin(); i!=Boundary.end(); i++)
		cout<<*i;
	
	vector<Edge> aux;
	set_difference(AllEdges.begin(),AllEdges.end(),Boundary.begin(),Boundary.end(),back_inserter(aux));
	cout<<"Lati interni:"<<endl;
	for(auto i=aux.begin(); i!=aux.end(); i++)
		cout<<*i;
}



void Grid::printedgesIndex(ofstream & ost1, ofstream & ost2, ofstream & ost3){
	set<unsigned int> s,ss;
	vector<unsigned int> diff;

	for (auto i : AllEdges) {s.insert(i.x()); s.insert(i.y());}
	ost1<<"All Edges:"<<endl;
	for (auto i : s) ost1<<i<<endl;

	for (auto i : Boundary) {ss.insert(i.x()); ss.insert(i.y());}
	ost2<<"Boundary Edges:"<<endl;
	for (auto i : ss) ost2<<i<<endl;

	//ostream_iterator<unsigned int> out(ost3,endl);
	//set_difference(s.begin(),s.end(),ss.begin(),ss.end(),out);
	set_difference(s.begin(),s.end(),ss.begin(),ss.end(),back_inserter(diff));
	ost3<<"Internal Edges:"<<endl;
	for (auto i : diff) ost3<<i<<endl;
	
}
*/
