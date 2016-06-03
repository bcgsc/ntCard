#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_

#include <fstream>

static inline bool
fileExists(std::string filename)
{
    std::ifstream in(filename.c_str());
	if (in.good()) {
		in.close();
		return true;
	} else {
		in.close();
		return false;
	}
}

#endif
