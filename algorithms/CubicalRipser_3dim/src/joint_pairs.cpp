/* joint_pairs.cpp

Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.

This file is part of CubicalRipser_3dim.

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
3 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>

#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "columns_to_reduce.h"
#include "simplex_coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"

using namespace std;

JointPairs2::JointPairs2(DenseCubicalGrids2* _dcg, ColumnsToReduce2* _ctr, vector<WritePairs2> &_wp, const bool _print){
	dcg = _dcg;
	ax = dcg -> ax;
	ay = dcg -> ay;
	az = dcg -> az;
	ctr = _ctr; // ctr is "0-dim"simplex list.
	ctr_moi = ctr -> max_of_index;
	n = ctr -> columns_to_reduce.size();
	print = _print;

	wp = &_wp;
	vtx = new Vertices2();

	for(int x = 1; x <= ax; ++x){
		for(int y = 1; y <= ay; ++y){
			for(int z = 1; z <= az; ++z){
				for(int type = 0; type < 3; ++type){
					int index = x | (y << 9) | (z << 18) | (type << 27);
					double birthday = dcg -> getBirthday2(index, 1);
					if(birthday < dcg -> threshold){
						dim1_simplex_list2.push_back(BirthdayIndex2(birthday, index, 1));
					}
				}
			}
		}
	}
	sort(dim1_simplex_list2.rbegin(), dim1_simplex_list2.rend(), BirthdayIndexComparator2());
}

void JointPairs2::joint_pairs_main2(){
	cubes_edges.reserve(2);
	UnionFind2 dset2(ctr_moi, dcg);
	ctr -> columns_to_reduce.clear();
	ctr -> dim = 1;
	double min_birth = dcg -> threshold;

	if(print == true){
		cout << "persistence intervals in dim " << 0 << ":" << endl;
	}
	
	for(auto e : dim1_simplex_list2){
		cubes_edges.clear();
		dcg -> GetSimplexVertices2(e.getIndex2(), 1, vtx);

		cubes_edges[0] = vtx -> vertex[0] -> getIndex2();
		cubes_edges[1] = vtx -> vertex[1] -> getIndex2();

		u = dset2.find2(cubes_edges[0]);
		v = dset2.find2(cubes_edges[1]);
			
		if(min_birth >= min(dset2.birthtime[u], dset2.birthtime[v])){
			min_birth = min(dset2.birthtime[u], dset2.birthtime[v]);
		}

		if(u != v){
			double birth = max(dset2.birthtime[u], dset2.birthtime[v]);
			double death = max(dset2.time_max[u], dset2.time_max[v]);

			if(birth == death){
				dset2.link2(u, v);
			} else {
				if(print == true){
					cout << "[" << birth << "," << death << ")" << endl;
				}
				wp -> push_back(WritePairs2(0, birth, death));
				dset2.link2(u, v);
			}
		} else { // If two values have same "parent", these are potential edges which make a 2-simplex.
			ctr -> columns_to_reduce.push_back(e);
		}
	}

	if(print == true){
		cout << "[" << min_birth << ", )" << endl;
	}

	wp -> push_back(WritePairs2(-1, min_birth, dcg -> threshold));
	sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator2());
}
