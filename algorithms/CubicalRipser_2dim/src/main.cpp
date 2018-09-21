/*
CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
2 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

#include "cubicalripser_2dim.h"

using namespace std;

int main(int argc, char** argv){
	const char* filename = nullptr;
	string output_filename = "answer_2dim.diagram"; //default name
	file_format format = DIPHA;
	calculation_method method = LINKFIND;
	double threshold = 99999;
	bool print = false;

	for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--threshold") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			threshold = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			string parameter = string(argv[++i]);
			if (parameter == "dipha") {
				format = DIPHA;
			} else if (parameter == "perseus") {
				format = PERSEUS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if(arg == "--method") {
			string parameter = string(argv[++i]);
			if (parameter == "link_find") {
				method = LINKFIND;
			} else if (parameter == "compute_pairs") {
				method = COMPUTEPAIRS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if (arg == "--output") {
			output_filename = string(argv[++i]);
		} else if(arg == "--print"){
			print = true;
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}
	
	ripser = CubicalRipser2D();
	ripser.ComputeBarcode(filename, output_filename, format, method, threshold, print);
	return 0;

}    