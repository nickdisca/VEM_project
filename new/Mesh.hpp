#ifndef _MESH_HPP_
#define _MESH_HPP_
#include <vector>
#include <string>
#include <memory>
#include <iosfwd>
/*
class Mesh; //pura dichiarazione, la definisco dopo

//serve puramente ad accedere agli elementi privati della classe
struct MeshHandler{
	MeshHandler(Mesh &mesh);
	//std::vector<Point> & pointList;
	//std::vector<Triangle> & elementList;
	//std::vector<Edge> & edgeList;
	//std::vector<Edge> & bEdgeList;
	//Mesh & m;
};


//class reader
class MeshReader{
public:
	MeshReader(bool verbose=false):M_verbose(verbose){}; //constructor
	~MeshReader(){}; //destructor
	// Reads a mesh
	int read(Mesh & mesh, std::string const & filename);
	// set verbosity
	void setverbose(bool set=true){M_verbose=set;};
protected:
	bool M_verbose;
};


//class mesh
class Mesh {
public:
	typedef std::size_t size_type;

	//standard methods + metodo che legge la mesh da file usando un reader
	Mesh()=default;
	Mesh(std::string const filename, MeshReader &);
	Mesh(const Mesh &)=default;
	Mesh(Mesh &&)=default;
	Mesh & operator =(Mesh const &)=default;
	Mesh & operator =(Mesh&&)=default;
	~Mesh()=default;

	friend struct MeshHandler;
	//! Number of points
	//size_type num_points()const {return M_pointList.size();}
	//! ith point
	//Point const & point(size_type i)const {return M_pointList[i];}
	//! Number of elements
	//size_type num_elements()const {return M_elementList.size();}
	//! ith element
	//Triangle const & element(size_type i)const {return M_elementList[i];}
	//! Ith element with specific name
	//Triangle const & triangle(size_type i)const {return M_elementList[i];}
	//! Number of edges
	//size_type num_edges()const {return M_edgeList.size();}
	//! ith Edge
	//Edge const & edge(size_type i)const {return M_edgeList[i];}
	//! Are edges stored?
	//bool has_Edges()const{return ! M_edgeList.empty();}
	//! Number of boundary edges
	//size_type num_bEdges()const {return M_bEdgeList.size();}
	//! ith edge
	//Edge const & bEdge(size_type i)const {return M_bEdgeList[i];}
	//! Are boundary edges stored
	//bool has_bEdges()const{return ! M_bEdgeList.empty();}
	//! Read a mesh from file using a reader
	//int readMesh(std::string const & file, MeshReader &);
	//! measure of the domain
	//double measure()const;
	//! Test mesh consistency
	//bool checkmesh()const;
private:
	//std::vector<Point> M_pointList;
	//std::vector<Triangle> M_elementList;
	//std::vector<Edge> M_edgeList;
	//std::vector<Edge> M_bEdgeList;
};
  
  //std::ostream & operator<<(std::ostream &, Mesh const&);




*/
#endif 
