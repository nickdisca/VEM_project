#ifndef _MESH_HPP_
#define _MESH_HPP_
#include <vector>
#include <string>
#include <memory>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <sstream>
#include <math.h>
#include "Geometry.hpp"
#include <Eigen/Sparse>

using MatrixType_S=Eigen::SparseMatrix<double>;

class Mesh; //pura dichiarazione, la definisco dopo

//serve puramente ad accedere agli elementi privati della classe
struct MeshHandler{
	MeshHandler(Mesh &mesh);
	std::vector<Point> & pointList;
	std::vector<Polygon> & elementList;
	std::vector<unsigned int> & boundary;
	unsigned int & num_edges;
	//std::vector<Edge> & edgeList;
	//std::vector<Edge> & bEdgeList;
	Mesh & m;
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
	Mesh(std::string const filename, MeshReader &, unsigned int kk);
	Mesh(const Mesh &)=default;
	Mesh(Mesh &&)=default;
	Mesh & operator =(Mesh const &)=default;
	Mesh & operator =(Mesh&&)=default;
	~Mesh()=default;

	friend struct MeshHandler;
	//! Number of points
	size_type num_points()const {return M_pointList.size();}
	//! ith point
	Point const & point(size_type i)const {return M_pointList[i];}
	//! Number of elements
	size_type num_elements()const {return M_elementList.size();}
	//! ith element
	Polygon const & element(size_type i)const {return M_elementList[i];}
	//! Ith element with specific name
	//Polygon const & triangle(size_type i)const {return M_elementList[i];}
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
	int readMesh(std::string const & file, MeshReader &, unsigned int kk);
	//! measure of the domain
	double area()const;
	double max_diam()const;
	//! Compute the connectivity matrix for the dof on edges
	void boundaryDOF();
	//! Test mesh consistency
	//bool checkmesh()const;
	//!Assemble global stiffness
	MatrixType_S GlobalStiffness(std::function<double (double, double)> mu, double mu_bar, bool constant_mu,
		std::function<double (double, double)> beta_x, std::function<double (double, double)> beta_y);
		
	MatrixType_S GlobalMass();
	MatrixType GlobalLoad(std::function<double(double,double)> f);

	std::vector<unsigned int> Dirichlet();
	MatrixType solve(std::function<double (double,double)> f, std::function<double (double,double)> g,
		std::function<double (double,double)> mu, double mu_bar, bool constant_mu, 
		std::function<double (double,double)> beta_x, std::function<double (double,double)> beta_y);

	//calcolo delle norme: da cambiare tutto con references
	MatrixType VEMConvert(std::function<double (double,double)> uex);
	
	double normInf(MatrixType & uex,MatrixType & u);
	void Allnorms(MatrixType & uex, MatrixType & u);
	
	//! output
	friend std::ostream & operator << (std::ostream &, const Mesh &);
private:
	std::vector<Point> M_pointList;
	std::vector<Polygon> M_elementList;
	std::vector<unsigned int> M_boundary;
	std::vector<Point> M_edgesDOF;
	unsigned int k;
	unsigned int M_num_edges;
	//std::vector<Edge> M_edgeList;
	//std::vector<Edge> M_bEdgeList;
};
  
  //std::ostream & operator<<(std::ostream &, Mesh const&);





#endif 
