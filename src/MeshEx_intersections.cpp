// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#include "MeshEx.h"
#include <list>

//#include <CGAL/Polygon_2.h>
//#include <CGAL/Triangulation_2.h>
//#include <CGAL/Triangulation_face_base_with_info_2.h>



//
//
//
bool MeshEx::findClosestIntersection( const math::Vec3d &position, const math::Vec3d &p1, const math::Vec3d &p2, math::Vec3d &intersection, math::Vec3d &normal, Triangle *&t )
{
	int count = 0;

	struct IntersectionHelper
	{
		IntersectionHelper( math::Vec3d intersection, Triangle *triangle, double sqDistance ) : m_intersection(intersection), m_triangle(triangle), m_sqDistance(sqDistance)
		{
		}

		// used for sorting
		bool operator<( const IntersectionHelper &rhs )
		{
			return m_sqDistance < rhs.m_sqDistance;
		}

		math::Vec3d m_intersection; // point of intersection
		Triangle       *m_triangle; // triangle in which the intersection point lies
		double         m_sqDistance; // distance to the near point
	};

	std::vector<IntersectionHelper> intersections;

	// test every triangle (TODO: do some spatial subdivision)
	for( std::vector<Triangle *>::iterator it = m_triangles.begin(); it != m_triangles.end(); ++it )
	{
		Triangle *tri = *it;

		Triangle_3 cgalTri( Point_3(tri->v0->position), Point_3(tri->v1->position), Point_3(tri->v2->position) );
		Segment_3 seg( Point_3(p1), Point_3(p2) );

		if( CGAL::do_intersect( cgalTri, seg ) )
		{
			CGAL::Point_3<K> tempIntersection;
			CGAL::Object o = CGAL::intersection(CGAL::Plane_3<K>( Point_3(tri->v0->position), Point_3(tri->v1->position), Point_3(tri->v2->position) ), seg);
			if( CGAL::assign(tempIntersection, o) )
			{
				double sd = CGAL::squared_distance( Point_3(position), tempIntersection );
				// intersection results in a point
				intersections.push_back( IntersectionHelper( math::Vec3d(tempIntersection.x(), tempIntersection.y(), tempIntersection.z()), tri, sd ) );
			}else
			{
				// segment!
				printf( "error : intersection during vertex movement is a line\n" );
			}
			//printf( "intersection\n" );
		}
	}



	// if we have found any intersection
	if( intersections.size() )
	{
		std::sort( intersections.begin(), intersections.end() );

		// find the closest one
		intersection = intersections[0].m_intersection;
		t = intersections[0].m_triangle;
		normal = t->normal;
		// and indicate success
		return true;
	}






	return false;
}