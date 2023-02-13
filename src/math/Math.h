// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#pragma once

#include <math.h>

#define MATH_PIf 3.14159265f
#define MATH_PI 3.14159265

#define MATH_2PIf 6.2831853f
#define MATH_INV_PIf 0.31830988618379067154f
#define MATH_INV_2PIf 0.15915494309189533577f
#define MATH_INV_4PIf 0.07957747154594766788f




#include "Vec2.h"
#include "Vec2Algo.h"
#include "Vec3.h"
#include "Vec3Algo.h"
#include "Matrix33.h"
#include "Matrix33Algo.h"
#include "Matrix44.h"
#include "Matrix44Algo.h"
#include "Color.h"
#include "RNG.h"
#include "Noise.h"

#include "PerlinNoise.h"

///
/// \brief vectors matrices and mathematical operations
///
namespace math
{
	inline int isnan( double f )
	{
		#if defined(_WINDOWS)
		return _isnan(f);
		#else
		return isnan(f);
		#endif
	}

	inline int isinf( double f )
	{
		#if defined(_WINDOWS)
		return !_finite((f));
		#else
		return isinf(f);
		#endif
	}

	inline double degToRad( float degree )
	{
		return (double) (( degree *  MATH_PIf)/180.0f );
	}

	inline double degToRad( double degree )
	{
		return (double) (( degree *  MATH_PI)/180.0 );
	}

	inline double sign( double f )
	{
		return f < 0.0f ? -1.0f : 1.0f;
	}


	//
	// color conversion functions
	//

	// RGBA

	#define REDSHIFT   0
	#define GREENSHIFT 8
	#define BLUESHIFT  16
	#define ALPHASHIFT 24

	inline unsigned long setColor( const unsigned int &r, const unsigned int &g, const unsigned int &b, const unsigned int &a )
	{
		return ((a << ALPHASHIFT) + (r << REDSHIFT) + (g << GREENSHIFT) + (b << BLUESHIFT));
	}

	//
	// Color related ops
	//

	inline Color operator+( const Color &lhs, const Color &rhs )
	{
	return Color( lhs.r+rhs.r, lhs.g+rhs.g, lhs.b+rhs.b, lhs.a+rhs.a );
	}
	inline Color operator-( const Color &lhs, const Color &rhs )
	{
	return Color( lhs.r-rhs.r, lhs.g-rhs.g, lhs.b-rhs.b, lhs.a-rhs.a );
	}
	inline Color operator*( const Color &lhs, const Color &rhs )
	{
	return Color( lhs.r*rhs.r, lhs.g*rhs.g, lhs.b*rhs.b, lhs.a*rhs.a );
	}
	inline Color operator/( const Color &lhs, const Color &rhs )
	{
	return Color( lhs.r/rhs.r, lhs.g/rhs.g, lhs.b/rhs.b, lhs.a/rhs.a );
	}

	inline Color operator+( const Color &lhs, const double &rhs )
	{
	return Color( lhs.r+rhs, lhs.g+rhs, lhs.b+rhs, lhs.a+rhs );
	}
	inline Color operator-( const Color &lhs, const double &rhs )
	{
	return Color( lhs.r-rhs, lhs.g-rhs, lhs.b-rhs, lhs.a-rhs );
	}
	inline Color operator*( const Color &lhs, const double &rhs )
	{
	return Color( lhs.r*rhs, lhs.g*rhs, lhs.b*rhs, lhs.a*rhs );
	}
	inline Color operator/( const Color &lhs, const double &rhs )
	{
	return Color( lhs.r/rhs, lhs.g/rhs, lhs.b/rhs, lhs.a/rhs );
	}

	inline Color operator+( const double &lhs, const Color &rhs )
	{
	return Color( (lhs + rhs.r), (lhs + rhs.g), (lhs + rhs.b), (lhs + rhs.a) );
	}
	inline Color operator-( const double &lhs, const Color &rhs )
	{
	return (rhs-lhs);
	}
	inline Color operator*( const double &lhs, const Color &rhs )
	{
	return (rhs*lhs);
	}
	inline Color operator/( const double &lhs, const Color &rhs )
	{
	return (rhs/lhs);
	}

	//
	// Inverts only the color channels of the given color, not the alpha channel
	//
	inline Color invert( const Color &col )
	{
		return Color( 1.0f - col.r, 1.0f - col.g, 1.0f - col.b, col.a );
	}




	double                                                   area( const Vec3d &p0, const Vec3d &p1, const Vec3d &p2 ); // computes area of an triangle

	double                                                                distance( const Vec3d &p0, const Vec3d &p1 ); // computes the euclidian distance between 2 points in space
	double                  distancePointPlane( const math::Vec3d &point, const Vec3d &normal, const double &distance ); // computes the distance of a point to a given plane


	math::Vec3d     projectPointOnPlane( const math::Vec3d &normal, const double &distance, const math::Vec3d &point ); // returns the projection of the given point on the normal and distance specified plane
	
	//
	// Misc mathematical utilities
	//

	//
	//
	//
	double mapValueTo0_1( const double &sourceRangeMin, const double &sourceRangeMax, const double &value );

	template<typename T>
	inline T max( T x, T y )
	{
		return x > y ? x : y;
	}

	template<typename T>
	inline T min( T x, T y )
	{
		return x < y ? x : y;
	}

}
