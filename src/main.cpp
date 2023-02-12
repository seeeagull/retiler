// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
#include <iostream>

#include "ObjIO.h"
#include "Mesh.h"
#include "MeshEx.h"
#include "RetileSurface.h"






int main(int argc, char ** argv)
{
	// load mesh
	dk::Mesh *m = dk::io::importFromObjFile( "../bunny.obj" );

	// generate winged edge datastructure
	MeshEx *mx = new MeshEx(m);

	delete m;

	// perform retiling
	int numVertices = 600; // manually specify the number of vertices (== level of detail)
	RetileSurface retiler;
	RetileSurface::ConstrainedVertexSet cvs; // this is a list of vertices which are fixed in their position (and dont follow repulsion)
	retiler.retile( mx, numVertices, cvs );

	// mx is now retiled and should have numVertices vertices etc.
	m = mx->getMesh();
	dk::io::exportToObjFile( m, "../simple_bunny.obj" );

	delete m;
	delete mx;

	return 0;
}
