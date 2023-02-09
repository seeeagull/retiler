// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#pragma once
#include "Vec3.h"
#include "Vec3Algo.h"
#include "Matrix33.h"
#include "BoundingBox3.h"
#include <limits>
#include <algorithm>
#include <cassert>



namespace math
{
	/// \brief a class which holds an origin and a target and is used by intersection routines
	template<typename T>
	class Ray3
	{
	public:
		Ray3();                                                                                   // constructor
		Ray3( const math::Vec3<T> &origin, const math::Vec3<T> &direction, const T &_tmin = std::numeric_limits<T>::min(), const T &_tmax = std::numeric_limits<T>::max());       // constructor

		math::Vec3<T>                                                    getPosition( T t )const; // returns origin+direction*t

		math::Vec3<T>                                                                     origin; // point in space where the ray originates from
		math::Vec3<T>                                                                  direction; // normalized direction of the ray
		T                                                                             tmin, tmax; // valid ray segment
	};

	

	// constructor
	template<typename T>
	Ray3<T>::Ray3() : tmin(0.0f), tmax(std::numeric_limits<T>::max())
	{
	}

	// constructor
	template<typename T>
	Ray3<T>::Ray3( const math::Vec3<T> &origin, const math::Vec3<T> &direction, const T &_tmin, const T &_tmax ) : origin(origin), direction(direction), tmin(_tmin), tmax(_tmax)
	{
	}


	// returns origin+direction*t
	template<typename T>
	math::Vec3<T> Ray3<T>::getPosition( T t )const
	{
		return origin + direction*t;
	}
	





	//
	// ray related functions
	//

	//bool RayHitSphere( const Vector &vSpherePosition, const double &fSphereRadius, CRay &oRay );
	//bool RayHitSphereValues( const Vector &vSpherePosition, const double &fSphereRadius, CRay &oRay, double &fHitNear, double &fHitFar, Vector *pHitPointNear, Vector *pHitPointFar );

	//bool RayHitPlane( const Vector &vPlaneNormal, const double &fD, CRay &oRay );
	//bool rayHitPlaneValues( const Vec3d &planeNormal, const doublee &planeDistance, Ray &raydoubleat &hitDistance, Vec3d *hitPoint );



	// computes a intersection between a ray and a plane
	//
	// returns true if an intersection occured, false otherwise
	template<typename T>
	bool intersectionRayPlane( const Ray3<T> &ray, const Vec3<T> &normal, const T &distance, Vec3<T> &hitPoint )
	{
		// project the ray direction onto the plane normal
		T temp = dot( ray.direction, normal );

		// if result is zero, then the direction is parallel to the plane -> no intersection
		if( !temp )
			return false;

		double hitDistance = -(dot( normal, ray.origin ) + distance) / temp;

		// the point must lie on the raysegment between origin and target to pass the test
		if( (hitDistance >= ray.tmax) || (hitDistance <= ray.tmin) )
			return false;

		hitPoint = ray.origin + hitDistance*ray.direction;

		return true;
	}


	// computes a intersection between a ray and another ray
	// algorithm based on Graphics Gems I page 304
	//
	// note: the 2 rays must be coplanar
	//
	// returns true if an intersection occured, false otherwise
	template<typename T>
	bool intersectionRayRay( const Ray3<T> &ray1, const Ray3<T> &ray2, Vec3<T> &hitPoint )
	{
		math::Vec3<T> cp = math::cross( ray1.direction, ray2.direction );
		T denom = cp.getSquaredLength();

		if( denom == (T)0.0 )
			// lines are parallel
			return false;

		// we need to compute s and t to test the line segments
		math::Vec3<T> c1 = ray2.origin - ray1.origin;

		T t = math::Matrix33<T>( c1.x, ray2.direction.x, cp.x, c1.y, ray2.direction.y, cp.y, c1.z, ray2.direction.z, cp.z ).getDeterminant() / denom;
		T s = math::Matrix33<T>( c1.x, ray1.direction.x, cp.x, c1.y, ray1.direction.y, cp.y, c1.z, ray1.direction.z, cp.z ).getDeterminant() / denom;


		// check line segments
		if( (t < ray1.tmin) || (s < ray2.tmin) || (t > ray1.tmax) || (s > ray2.tmax) )
			return false;

		// compute intersection point
		hitPoint = ray2.origin + ray2.direction*s;

		// check for coplanarity
		//if( hitPoint != (ray1.originrigin + ray1.directionirection*t) )
		//	return false;

		// done
		return true;
	}


	// sphere is assumed to be at origin, r is radius
	template<typename T>
	bool intersectionRaySphere( const Ray3<T> &ray, T radius, T &t0, T &t1 )
	{
		math::Vec3<T> vec = ray.origin;
		T a = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y + ray.direction.z * ray.direction.z;
		T b = 2 * (ray.direction.x * vec.x +  ray.direction.y * vec.y + ray.direction.z * vec.z);
		T c = vec.x * vec.x + vec.y * vec.y + vec.z * vec.z - radius * radius;

		if (b == 0)
		{
			// handle special case where the the two vector ray.directionir and V are perpendicular
			// with V = ray.originrig - sphere.centre
			if (a == 0) return false;
			t0 = 0; t1 = sqrt(-c/a);
		}else
		{
			T discr = b * b - 4 * a * c;
			if (discr < 0) return false;
			T q = (b < 0) ? (T)-0.5 * (b - sqrt(discr)) : (T)-0.5 * (b + sqrt(discr));
			t0 = q / a;
			t1 = c / q;
		}


		if (t0 > t1) std::swap(t0, t1);
		return true;
	}

	template<typename T>
	bool intersectionRayBox( const Ray3<T> &ray, const BoundingBox3<T> &box, T &hitt0, T &hitt1 )
	{
		T t0 = ray.tmin, t1 = ray.tmax;
		for (int i = 0; i < 3; ++i)
		{
			// Update interval for _i_th bounding box slab
			T invRayDir = 1.f / ray.direction[i];
			T tNear = (box.minPoint[i] - ray.origin[i]) * invRayDir;
			T tFar  = (box.maxPoint[i] - ray.origin[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$s
			if (tNear > tFar) std::swap(tNear, tFar);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar  < t1 ? tFar  : t1;
			if (t0 > t1) return false;
		}
		hitt0 = t0;
		hitt1 = t1;
		return true;
	}




	//bool RayHitTriangle( const Vector &vPoint1, const Vector &vPoint2, const Vector &vPoint3, CRay &oRay );
	//bool RayHitTriangleValues( const Vector &vPoint1, const Vector &vPoint2, const Vector &vPoint3, CRay &oRay, double &fHit, Vector *pHitPoint,doublet *pUdoubleat *pV );

	typedef Ray3<float> Ray3f;
	typedef Ray3<double> Ray3d;
}
