// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This class performs the algorithm which is proposed in the paper
"Retiling Polygonal Surfaces" from G.Turk. The Algorithm uses
the CGAL for computational geometry stuff.

----------------------------------------------------------------------*/
#pragma once

#include "MeshEx.h"
#include <stdarg.h>




///
/// \brief performs retiling like proposed in Turk'92
///
class RetileSurface
{
public:
	typedef std::vector<std::pair<math::Vec3d, MeshEx::Triangle *> > ConstrainedVertexSet;
	RetileSurface();                                                                                         // constructor
	~RetileSurface();                                                                                        // destructor
	void                             retile( MeshEx *_mesh, size_t vertexCount, ConstrainedVertexSet &cvs ); // performs the algorithm on the given mesh with a set of contrained vertices





private:
	//
	// The vertex structure is used to create the positions of the new vertices and move them into
	// the right place
	//
	struct Vertex
	{
		Vertex( double x, double y, double z )
		{
			position = math::Vec3d( x, y, z );
			force = math::Vec3d( 0.0f, 0.0f, 0.0f );
			fixed = false;

			tag = false;
		}
		math::Vec3d position;

		MeshEx::Triangle                  *t; // reference to the triangle this vertex lies on
		math::Vec3d                    force; // superposition of all repulsing forces from neighbouring points
		bool                           fixed; // if this member is set to true, then the vertex wont be moved during relaxation

		std::vector<math::Vec3d>        path;



		bool tag;
		size_t index;
	};

	void                                                                 generateVertices( size_t number ); // generates a specified number of vertices which will be randomly distributed (not evenly) over the mesh
	void                                                                                   clearVertices(); // deletes and removes all vertices

	void                                                      performIteration( double radius, double damp ); // this method moves the vertices over the mesh surface, so that they are evenly distributed
	math::Vec3d mapPointToPlane( const math::Vec3d &q, MeshEx::Triangle *t, MeshEx::Edge *hingeEdge = 0  ); // maps a neighbour of a vertex onto a plane - if the neighbour lies on a neighbouring tri, then it is rotated around the shared edge
	void                                  moveVertexOnMesh( Vertex *vertex, const math::Vec3d &direction ); // moves the vertex over the meshsurface into the given direction

	void                                                                             doMutualTesselation(); // after the new vertices have been created and moved, they are triangulated into the original mesh

	void           writeVertices( const std::vector<Vertex *> &vertices, const char *filenameFormat, ... ); // utility method for debug and visualization-dumps out the positions of the vertex candidates into a binary file

	MeshEx                                                                                           *mesh; // mesh on which the algorithm will perform

	std::vector<Vertex *>                                                                         vertices; // the new vertices which will be used

	double                                                                                         m_radius; // radius which is used to get neighbouring vertices (will be computed)
	double                                                                                           m_damp; // damping factor of the repulsion forces
	size_t                                                                                m_iterationCount; // number of iterations for moving the vertices into their place
};




