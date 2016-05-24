#include "grid.hpp"
#include <string>

using namespace Geometry;
using namespace std;

Grid::Grid(ifstream & file){
string str,aux;
unsigned int Npoly=0,Nvert=0,tipo=0;
set<Edge> insieme, interni;
Edge e;
double i1,i2;
vector<unsigned int> lati;
pair<set<Edge>::iterator,bool> check;

//read the first line
getline(file,str);
stringstream ss(str);
//getline(ss,aux,' '); Nvert=stoi(aux); getline(ss,aux,' '); Npoly=stoi(aux);
ss>>Nvert>>Npoly;
//cout<<Nvert<<"   "<<Npoly<<endl;

vect.resize(Nvert);

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
	lati.clear();
	while (getline(ss,aux,' ')){
		//vvv.emplace_back(vect[stoi(aux)]);
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
		case(0): {Triangle t(lati, &vect); abspol.emplace_back(make_shared<Triangle> (t)); t.showMe(); cout<<"Area"<<t.area()<<endl;} break;
		case(1): {Square s(lati, &vect); abspol.emplace_back(make_shared<Square> (s)); s.showMe(); cout<<"Area"<<s.area()<<endl;} break;
		default: {Polygon p(lati, &vect); abspol.emplace_back(make_shared<Polygon> (p)); p.showMe(); cout<<"Area"<<p.area()<<endl;}
	}
}

//copy all edges in the vector
AllEdges.resize(insieme.size());
copy(insieme.begin(),insieme.end(),AllEdges.begin());

//create the vector containing Boundary edges
Boundary.resize(insieme.size()-interni.size());
set_difference(insieme.begin(), insieme.end(),interni.begin(),interni.end(),Boundary.begin());
//set_difference(insieme.begin(), insieme.end(),interni.begin(),interni.end(),back_inserter(Boundary));
}


double Grid::area(){
	unsigned int Npoly=abspol.size();
	double area{0.0};
	for (unsigned int i=0; i<Npoly; i++){
		area+=abspol[i]->area();
	}
return area;
}



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