// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#include "Math.h"




namespace math
{

	//
	// computes area of an triangle
	//
	double area( const Vec3d &p0, const Vec3d &p1, const Vec3d &p2 )
	{
		double la = (p1 - p0).getLength(); // compute lengths of the triangle sides
		double lb = (p2 - p1).getLength();
		double lc = (p2 - p0).getLength();
		double s = 0.5f*( la+lb+lc ); // compute the semiperimeter
		return sqrt( s*(s-la)*(s-lb)*(s-lc) ); // compute the area
	}


	//
	// computes the distance of a point to a plane
	//
	inline double distancePointPlane( const math::Vec3d &point, const Vec3d &normal, const double &distance )
	{
		return dotProduct( normal, point ) + distance;
	}

	//
	// computes the euclidian distance between 2 points in space
	//
	double distance( const Vec3d &p0, const Vec3d &p1 )
	{
		return (p1-p0).getLength();
	}


	//
	// returns the projection of the given point on the normal and distance specified plane
	//
	math::Vec3d projectPointOnPlane( const math::Vec3d &normal, const double &distance, const math::Vec3d &point )
	{
		return point - distancePointPlane( point, normal, distance )*normal;
	}


	//
	//
	//
	double mapValueTo0_1( const double &sourceRangeMin, const double &sourceRangeMax, const double &value )
	{
		return (value-sourceRangeMin) / (sourceRangeMax - sourceRangeMin);
	}

}
