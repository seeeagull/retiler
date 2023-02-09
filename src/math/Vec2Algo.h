// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#pragma once

#include "Vec2.h"
#include "Vec2Algo.h"


namespace math
{


	template<typename T>
	inline Vec2<T> operator-( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x-rhs.x, lhs.y-rhs.y );
	}

	template<typename T>
	inline Vec2<T> operator+( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x+rhs.x, lhs.y+rhs.y );
	}

	template<typename T>
	inline Vec2<T> operator*( const Vec2<T> &lhs, const double &rhs )
	{
		return Vec2<T>( lhs.x*rhs, lhs.y*rhs );
	}

	template<typename T>
	inline Vec2<T> operator*( const T &lhs, const Vec2<T> &rhs )
	{
		return (rhs*lhs);
	}

	template<typename T>
	inline Vec2<T> operator*( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x*rhs.x, lhs.y*rhs.y );
	}

	template<typename T>
	inline Vec2<T> operator-( const Vec2<T> &rhs )
	{
		return Vec2<T>( -rhs.x, -rhs.y );
	}

	template<typename T>
	inline double crossProduct( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return lhs.x*rhs.y - lhs.y*rhs.x;
	}

	template<typename T>
	inline void crossProduct( double &result, const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		result = lhs.x*rhs.y - lhs.y*rhs.x;
	}



}
