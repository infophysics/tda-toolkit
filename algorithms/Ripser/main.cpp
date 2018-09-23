


#include "ripser.h"
#include <string>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>


int main(int argc, char** argv) {

	const char* filename = nullptr;
	long dim_max = 1;
	float threshold = 0;
	std::string format;
	float ratio = 1;
	
	Ripser ripser;

	

#ifdef USE_COEFFICIENTS
	coefficient_t modulus = 2;
#else
	const coefficient_t modulus = 2;
#endif
	for (index_t i = 1; i < argc; ++i) {
			const std::string arg(argv[i]);
			if (arg == "--help") {
				ripser.print_usage_and_exit(0);
			} else if (arg == "--dim") {
				std::string parameter = std::string(argv[++i]);
				size_t next_pos;
				dim_max = std::stol(parameter, &next_pos);
				if (next_pos != parameter.size()) ripser.print_usage_and_exit(-1);
			} else if (arg == "--threshold") {
				std::string parameter = std::string(argv[++i]);
				size_t next_pos;
				threshold = std::stof(parameter, &next_pos);
				if (next_pos != parameter.size()) ripser.print_usage_and_exit(-1);
			} else if (arg == "--ratio") {
				std::string parameter = std::string(argv[++i]);
				size_t next_pos;
				ratio = std::stof(parameter, &next_pos);
				if (next_pos != parameter.size()) ripser.print_usage_and_exit(-1);
			} else if (arg == "--format") {
				std::string parameter = std::string(argv[++i]);
				if (parameter == "lower-distance")
					format = parameter;
				else if (parameter == "upper-distance")
					format = parameter;
				else if (parameter == "distance")
					format = parameter;
				else if (parameter == "point-cloud")
					format = parameter;
				else if (parameter == "dipha")
					format = parameter;
				else if (parameter == "ripser")
					format = parameter;
				else
					ripser.print_usage_and_exit(-1);
	#ifdef USE_COEFFICIENTS
			} else if (arg == "--modulus") {
				std::string parameter = std::string(argv[++i]);
				size_t next_pos;
				modulus = std::stol(parameter, &next_pos);
				if (next_pos != parameter.size() || !is_prime(modulus)) ripser.print_usage_and_exit(-1);
	#endif
			} else {
				if (filename) { ripser.print_usage_and_exit(-1); }
				filename = argv[i];
				
			}
	}
	ripser.ComputeBarcode(filename, dim_max, threshold, ratio, format, modulus);
	return 0;
}
