// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#pragma once



namespace math
{
	///
	/// \brief similar to the math::Vec3d but specialized to usage as color
	///
	class Color
	{
	public:
		Color();
		Color( const double &r, const double &g, const double &b, const double &a = 1.0f );
		~Color();

		// standard colors
		static Color                                                              White();
		static Color                                                              Black();
		static Color                                                               Blue();
		static Color                                                             Yellow();
		static Color                                                              Green();
		static Color                                                                Red();
		static Color From255( const unsigned char &r, const unsigned char &g, const unsigned char &b, const unsigned char &a = 255 );


        void set( const double &r, const double &g, const double &b, const double &a = 1.0f ); 

		void                                                                clamp( void ); ///< clamp the component values into the range of [0,1]
        void                                                               invert( void ); ///< invert the color -> each component is 1.0f - value
		unsigned long                                                   makeDWORD( void ); ///< returns a RGBA(4x8bit) representation of the color

		// operators
		bool                                               operator==( const Color &rhs );
        bool                                               operator!=( const Color &rhs );
        bool                                               operator+=( const Color &rhs );
        bool                                               operator-=( const Color &rhs );
        bool                                               operator*=( const Color &rhs );

        bool                                               operator+=( const double &rhs );
        bool                                               operator-=( const double &rhs );
        bool                                               operator*=( const double &rhs );
        bool                                               operator/=( const double &rhs );

        const                                            double& operator[]( int i ) const; ///< returns ith component
        double&                                                        operator[]( int i ); ///< returns ith component

		union
		{
			struct
			{
				double	r, g, b, a;
			};
			double v[4];
		};
	};
}
