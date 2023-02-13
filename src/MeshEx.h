// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
#pragma once
#include <map>
#include <vector>


#include "math/Math.h"

// CGAL -------------------------
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Filtered_kernel< CGAL::Simple_cartesian<double> > K;

typedef CGAL::Segment_3<K>                              Segment_3;
typedef CGAL::Triangle_3<K>                            Triangle_3;
CGAL::Point_3<K>                Point_3( const math::Vec3d &vec );

#include "Mesh.h"

using namespace dk;


///
/// \brief winged edge data structure
///
class MeshEx
{
public:
	struct Triangle;

	///
	/// \brief correspondeces to the vertices of the input mesh
	///
	struct Vertex
	{
		Vertex();                                       // constructor
		void                                               registerTriangle( Triangle *t );  // adds the given element to the elementring list
		void                                             unRegisterTriangle( Triangle *t );  // remvoes the given element to the elementring list
		bool                                                                  isDesolate();  // returns true if no element references this node
		void                                                               computeNormal();  // computes normal from the surrounding tris - works only if the surrounding tris have a valid normal

		math::Vec3d                                                               position;
		math::Vec3d                                                                 normal;
		math::Color                                                                  color;

		std::vector<Triangle *>                                               triangleRing;  // references to all elements which index into this node
	};

	///
	/// \brief We use the winged edge datastructure to ease mesh edit operations
	///
	struct Edge
	{
		Edge( Vertex *_v1, Vertex *_v2 );                                                                            // constructor

		void                                                                       registerTriangle( Triangle *t );  // assigns the given element to the left or right wing of the edge (dependand on which wing is 0)
		void                                                                     unregisterTriangle( Triangle *t );  // sets either the left or right reference to zero when it references the given element
		Vertex                                                                        *getOtherVertex( Vertex *v );  // convinience function to get the other node of the two nodes which are referenced by the edge
		Triangle                                                                  *getOtherTriangle( Triangle *t );  // returns the opposite triangle
		bool                                                                                          isDesolate();  // returns true if no elements references this edge
		bool                                                                                 contains( Vertex *v );  // returns true if the edge contains the given vertices
		bool                                                          contains( Vertex *vertex1, Vertex *vertex2 );  // returns true if the edge contains the given vertices
		bool                                                                               contains( Triangle *t );  // returns true if the edge contains the given triangles
		bool                                                                                      isBoundaryEdge();  // returns true if the edge has only one element reference instead of two

		Vertex                                                                                                 *v1;  // Edge is made up of 2 nodes
		Vertex                                                                                                 *v2;
		Triangle                                                                                             *left;  // left neighbour - if it is 0, the edge is a borderedge
		Triangle                                                                                            *right;  // right neighbour - if it is 0, the edge is a borderedge

		bool tag;
	};

	///
	/// \brief correspondences to the triangles of the input mesh
	///
	struct Triangle
	{
		union
		{
			struct
			{
				Vertex *v0, *v1, *v2;
			};
			Vertex *v[3];
		};
		union
		{
			struct
			{
				Edge *e1, *e2, *e3;
			};
			Edge *e[3];
		};

		Triangle();                                                                                         // constructor
		void                                                                           computeProperties(); // computes normal, area, etc...
		Vertex                                         *getOtherVertex( Vertex *vertex1, Vertex *vertex2 ); // returns the third node which neither equals n1 nor n2
		Vertex                                                               *getOtherVertex( Edge *edge ); // returns the third node which neither equals n1 nor n2
		bool                                                              contains( Vertex *vertex ) const; // returns true if one of the triangles vertices equals the given vertex
		size_t                                                      getVertexsLocalIndex( Vertex *vertex ); // returns the local index of the node within the element
		Edge                                                                 *getSharedEdge( Triangle *t ); // returns the edge which is shared with the given triangle or 0 if they dont share an edge
		Edge                                                  *getEdge( Vertex *vertex1, Vertex *vertex2 ); // returns the edge of the triangle which contains the 2 given vertices
		Edge                                                                    *getOtherEdge( Vertex *v ); // returns the one edge of the triangle which does not contain the vertex v

		void                                                                               computeNormal(); // (re) computes normal

		math::Vec3d                                                                                 center; // barycentric center position of the triangle
		math::Vec3d                                                                                 normal;
		math::Vec3d                                                                                      u; // these 2 vectors (u and v) define the local 2dimensional coordinate
		math::Vec3d                                                                                    v_v; // frame of the element (although they are 3d vectors)

		math::Matrix33f                                                                               beta; // barycentric coordinate matrix -> base vectors are the position of the 3 nodes of the element in the local 2d cooridnate frame
																											// so beta transforms local coordinates into barycentric coordinates

		double                                                                                         area; // the area of the element

		bool tag;
		bool tag2;
	};


	MeshEx( Mesh *mesh );                                                                             ///< constructor which takes a mesh as input
	~MeshEx();                                                                                        ///< destructor
	Mesh                                                                                  *getMesh(); ///< creates a mesh which represents the MeshEx


	Vertex                                              *createVertex( const math::Vec3d &position ); ///< creates a Vertex and adds it to the node list with the given world position
	void                                                              removeVertex( Vertex *vertex ); ///< removes given node from the node list
	Triangle     *createTriangle( const size_t &index0, const size_t &index1, const size_t &index2 ); ///< creates a Triangle and adds it to the element list with the given node indices
	Triangle                                   *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2 ); ///< creates a Triangle and adds it to the triangle list from given vertices
	Triangle     *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2, Edge *e0, Edge *e1, Edge *e2 ); ///< creates a Triangle and adds it to the triangle list from given edges
	Triangle       *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2, std::vector<Edge *> &edges ); ///< creates a Triangle and adds it to the triangle list -> edges are looked up within the given edges (for speedup)
	Edge                                                       *createEdge( Vertex *v1, Vertex *v2 ); ///< creates a new edge from 2 given m_vertices
	Edge                                                         *findEdge( Vertex *n1, Vertex *n2 ); ///< looks for an edge between the given nodes - returns 0 if it could not be found
	Vertex *                               splitEdge( Edge *edge, const math::Vec3d &splitPosition ); ///< splits the given edge into 2 and returns the splitnode which was created to seperate the edge
	void                                                                    removeEdge( Edge *edge ); ///< removes edge from the list of edges
	void    removeTriangle( Triangle *element, bool removeVertices = true, bool removeEdges = true ); ///< removes element from the list of elements
	void                  reTriangulate( Triangle *t, const std::vector<math::Vec3d> &pointsWithin ); ///< retriangulates the given triangle taking the given points into acount (which have to be coplanar with the triangle)
	bool                                      removeVertexAndReTriangulateNeighbourhood( Vertex *v ); ///< removes the given vertex and retriangulate the hole which would be made
	void           createTrianglesFromEdges( std::vector<Edge *> &edges, const math::Vec3d &normal ); ///< this method takes a list of already created edges and creates triangles for filling triangulated areas

	void retriangulateHole( std::vector<MeshEx::Edge *> &boundaryEdges, std::map<MeshEx::Vertex *, math::Vec2d> &boundaryVertexProjections, std::vector<std::pair<math::Vec3d, math::Vec2d> > &interiorPoints );
	void                                                                        detectAndFillHoles();

	void                                                                                     clear();


	std::vector<Vertex*>                                                                  m_vertices;
	std::vector<Triangle*>                                                               m_triangles;
	std::vector<Edge*>                                                                       m_edges;

	bool     stop;
	math::Vec3d temp;
};
