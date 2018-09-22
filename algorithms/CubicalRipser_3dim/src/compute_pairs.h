/* compute_pairs.h

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


#include <vector>
#include <unordered_map>

using namespace std;

template <class Key, class T> class hash_map2 : public std::unordered_map<Key, T> {};

class ComputePairs2
{
public:
	DenseCubicalGrids2* dcg;
	ColumnsToReduce2* ctr;
	hash_map2<int, int> pivot_column_index;
	int ax, ay, az;
	int dim;
	vector<WritePairs2> *wp;
	bool print;

	ComputePairs2(DenseCubicalGrids2* _dcg, ColumnsToReduce2* _ctr, vector<WritePairs2> &_wp, const bool _print);

	void compute_pairs_main2();

	void outputPP2(int _dim, double _birth, double _death);

	BirthdayIndex2 pop_pivot2(priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndexComparator2>&
		column);

	BirthdayIndex2 get_pivot2(priority_queue<BirthdayIndex2, vector<BirthdayIndex2>, BirthdayIndexComparator2>&
		column);

	void assemble_columns_to_reduce2();
};