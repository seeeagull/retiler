// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This header offers some very usefull functions which you can use
for manipulating and working with std::strings.

----------------------------------------------------------------------*/
#include "StringManip.h"

namespace dk
{
	namespace util
	{
		inline std::string trim_right (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			std::string::size_type i (d.find_last_not_of (t));
			if (i == std::string::npos)
				return "";
			else
				return d.erase (d.find_last_not_of (t) + 1) ; 
		}  // end of trim_right

		inline std::string trim_left (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			return d.erase (0, s.find_first_not_of (t)) ; 
		}  // end of trim_left

		inline std::string trim (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			return trim_left (trim_right (d, t), t) ; 
		}  // end of trim


		// split a line into the first word, and rest-of-the-line
		std::string getWord (std::string & s, const std::string delim,const bool trim_spaces)
		{
			// find delimiter  
			std::string::size_type i (s.find (delim));

			// split into before and after delimiter
			std::string w (s.substr (0, i));

			// if no delimiter, remainder is empty
			if (i == std::string::npos)
				s.erase ();
			else
				// erase up to the delimiter
				s.erase (0, i + delim.size ());

			// trim spaces if required
			if (trim_spaces)
			{
				w = trim (w);
				s = trim (s);
			}

			// return first word in line
			return w;
		} // end of getWord	

		// To be symmetric, we assume an empty string (after trimming spaces)
		// will give an empty vector.
		// However, a non-empty string (with no delimiter) will give one item
		// After that, you get an item per delimiter, plus 1.
		// eg.  ""      => empty
		//      "a"     => 1 item
		//      "a,b"   => 2 items
		//      "a,b,"  => 3 items (last one empty)
		void splitString( const std::string s, std::vector<std::string> & v, const std::string delim, const bool trim_spaces)
		{
			// start with initial string, trimmed of leading/trailing spaces if required
			std::string s1 (trim_spaces ? trim (s) : s);

			v.clear (); // ensure vector empty
			
			// no string? no elements
			if (s1.empty ())
				return;

			// add to vector while we have a delimiter
			while (!s1.empty () && s1.find (delim) != std::string::npos)
				v.push_back (getWord (s1, delim, trim_spaces));

			// add final element
			v.push_back (s1);
		} // end of splitString 


	} // namespace util
} // namespace dk