#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "Cells/All.h"
#include "Complexes/All.h"
#include "Algos/All.h"
#include "Perseus.h"

int main(int argc, char* argv[])
{
	
	if (argc < 3 || argc > 5)
	{
		cout<<"\nIncorrect number of arguments!";
		cout<<"\nUsage: <Program Name> <Input Type> <Input Filename> (optional)<Output Filename> (optional)<Engine>\n";
		cout<<"Where: \n  1. The acceptable <Input Type> is SimTop, CubTop, Rips, etc... \n";
		cout<<"  2. If <Output Filename> is unspecified, \"output\" will be used with suffixes.\n";
		cout<<"  3. The legal <engine> values are R, C and A for reduction, coreduction and alternation.";
		return -1;

	}

	// obtain string representations of input options:
	string input_type = argv[1];
	string infile = argv[2];
	string outfile = (argc >= 4)? argv[3] : "output";
	string eng = (argc >= 5) ? argv[4] : "a";

	//cout<<"\neng is "<<eng; cin.get();
	pers = Perseus();
	pers.ComputeBarcode(input_type, infile, outfile, eng);
}