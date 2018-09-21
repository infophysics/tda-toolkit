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

#include "dense_cubical_grids.h"
#include "birthday_index.h"
#include "columns_to_reduce.h"
#include "simplex_coboundary_enumerator.h"
#include "write_pairs.h"
#include "union_find.h"
#include "compute_pairs.h"

using namespace std;

enum calculation_method { LINKFIND, COMPUTEPAIRS};

class CubicalRipser2D{
	
	public:
		CubicalRipser2D();
		virtual ~CubicalRipser2D();
		void print_usage_and_exit(int exit_code);
		void ComputeBarcode(const char* filename, string output_filename, string format, string method, double threshold, bool print);	
};


