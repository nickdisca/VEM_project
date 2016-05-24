#ifndef HH_EDGE_HH
#define HH_EDGE_HH
#include "Polygon.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iterator>

using namespace Geometry;
using namespace std;

class Edge {
public:
	//Edge ():a(0),b(0){};
	Edge()=default;
    Edge (const Edge & )=default;
    Edge & operator= (const Edge & )= default;
    Edge(unsigned int aa, unsigned int bb): a(aa),b(bb){};

    void set(unsigned int const &aa, unsigned int const &bb){
      a=aa; b=bb;};
    unsigned int x(){return a;};
    unsigned int y(){return b;};
    void print(){cout<<a<<" "<<b<<endl;};
    friend bool operator <(Edge const &f, Edge const &s)
    {
    	auto mm=minmax(f.a,f.b);
    	auto mm2=minmax(s.a,s.b);
    	if (mm.first==mm2.first) return mm.second<mm2.second;
    	return mm.first<mm2.first;
    };
    friend ostream & operator << (ostream & ost, Edge const & e){
    	ost<<e.a<<" "<<e.b<<endl;
    	return ost;
    };



private:
	unsigned int a;
	unsigned int b;
};



#endif