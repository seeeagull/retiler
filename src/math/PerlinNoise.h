// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
//
//
// PerlinNoise generator (1d,2d,3d,4d)
//
//
#pragma once


namespace math
{
	//
	// this class is a simple utilityclass implementing perlin noise
	//
	class PerlinNoise
	{
	public:

		PerlinNoise();     // constructor - initializes perlin noise parameters here

		double                    perlinNoise_2D( double s, double t ); // 2D perlin noise function
		double          perlinNoise_3D( double u, double v, double w  ); // 3D perlin noise function
		double perlinNoise_4D( double u, double v, double w, double x  ); // 4D perlin noise function

		double                                        getAmplitude( void );
		void                              setAmplitude( double amplitude );
		double                                   getAmplitudeRatio( void );
		void                    setAmplitudeRatio( double amplitudeRatio );
		double                                        getFrequency( void );
		void                              setFrequency( double frequency );
		double                                   getFrequencyRatio( void );
		void                    setFrequencyRatio( double frequencyRatio );
		int                                              getDepth( void );
		void                                        setDepth( int depth );
		bool                                        getInflection( void );
		void                             setInflection( bool inflection );

	private:
		// perlin noise parameters
		double                                                 m_amplitude;
		double                                            m_amplitudeRatio; // will control the influence of the higherfrequencies
		double                                                 m_frequency;
		double                                            m_frequencyRatio;

		int                                                     m_octaves; // number of passes - the more passes, the more higher frequencies
		bool                                                 m_inflection;

		static unsigned char                    g_permutationTable[ 256 ];
		static double                               g_gradientTable[ 768 ];
		static bool                                     g_tablesInitiated;

		static void                         initGradientTable( int seed );

        double                       interpolatedNoise( double s, double t ); // this is the noise interpolation function for 2dimensional perlin noise
        double              interpolatedNoise( double u, double v, double w ); // this is the noise interpolation function for 3dimensional perlin noise
        double interpolatedNoise( double _u, double _v, double _w, double _s ); // this is the noise interpolation function for 4dimensional perlin noise

		double               interpolatedGradientNoise( double u, double v );
		double      interpolatedGradientNoise( double u, double v, double w );
		double               interpolatedGradientNoise( double _u, double _v,
			                                         double _w, double _s );

		double                                              noise( int x ); // this is a noise function which will return a pseudorandom number based on the value x
		double                                       noise( int x, int y );
		double                                noise( int x, int y, int z );
		double                         noise( int x, int y, int z, int w );

		double            gradientNoise(int x, int y, double fx, double fy );
		double                           gradientNoise(int x, int y, int z,
			                                double fx, double fy, double fz);
		double                    gradientNoise(int x, int y, int z, int w,
			                     double fx, double fy, double fz, double fw );

		double            interpolateLinear( double v1, double v2, double t ); // simple linear interpolation between two values dependent on a interpolation value
        double         interpolateEaseCurve( double v1, double v2, double t ); // will interpolate between two values in a more harmonic and smooth manner
	};
}
