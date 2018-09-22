/* simplex_coboundary_enumerator.h

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


class SimplexCoboundaryEnumerator2
{
public:
	BirthdayIndex2 simplex;
	DenseCubicalGrids2* dcg;
	Vertices2* vtx;
	double birthtime;
	int ax, ay, az;
	int cx, cy, cz;
	int count;
	BirthdayIndex2 nextCoface;
	double threshold;

	SimplexCoboundaryEnumerator2();

	void setSimplexCoboundaryEnumerator2(BirthdayIndex2 _s, DenseCubicalGrids2* _dcg); 

	bool hasNextCoface2(); 

	BirthdayIndex2 getNextCoface2();
};