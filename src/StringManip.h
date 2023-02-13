// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This header offers some very usefull functions which you can use
for manipulating and working with std::strings.

----------------------------------------------------------------------*/
#pragma once
#include <string>
#include <vector>

namespace dk
{
	namespace util
	{
		inline std::string trim_right (const std::string & s, const std::string & t = " \t\r\n");

		inline std::string trim_left (const std::string & s, const std::string & t = " \t\r\n");

		inline std::string trim (const std::string & s, const std::string & t = " \t\r\n");


		/// split a line into the first word, and rest-of-the-line
		std::string getWord (std::string & s, const std::string delim = " ",const bool trim_spaces = true);

		/// split string through some tokens
		void splitString( const std::string s, std::vector<std::string> & v, const std::string delim = " ", const bool trim_spaces = true);

	} // namespace util
} // namespace dk










