#pragma once
/*
 * Pers.cpp
 *
 * Contains main() function for Persistent Homology
 * Using Discrete Morse Theory
 */
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <unistd.h>

#include "Cells/All.h"
#include "Complexes/All.h"
#include "Algos/All.h"

// class from which the birth times arise, usually = "int"
typedef int BC;
// class defining the ring over which we perform chain algebra, usually = "int"
typedef int CC;

using namespace std;

class Perseus{
	
	public:
		Perseus();
		virtual ~Perseus();
		// silly string functions for dealing with input arguments, change case
		inline void lowerCase(std::string& s)
		{std::transform(s.begin(), s.end(), s.begin(), ::tolower);}

		inline void upperCase(std::string& s)
		{std::transform(s.begin(), s.end(), s.begin(), ::toupper);}

		// wrapper declaration for converting input file data to cell complexes;
		// actual functions follow main():

		// sparse cubical complex from top cube info
		bool buildCToplexFromFile (ifstream&, MComplex<CC,BC>&, CToplex<CC,int,BC>&);
		// dense cubical complex from top cube info
		bool buildDenseCToplexFromFile (ifstream&, MComplex<CC,BC>&, DenseCToplex<CC,BC>&);
		// uniform triangulation from top simplex info
		bool buildSToplexFromFile (ifstream&, MComplex<CC,BC>&, SToplex<CC,double,BC>&, bool, bool);
		// non-uniform triangulation from top simplex info
		bool buildNMFSToplexFromFile (ifstream&, MComplex<CC,BC>&, SToplex<CC,double,BC>&, bool, bool);
		// rips complex from points with non-uniform birth times
		bool buildRIPSFromFile (ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&);
		// typical rips complex, growing balls around points
		bool buildBRIPSFromFile (ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&,bool,bool,bool);
		// rips complex from point-point distance matrix
		bool buildBRIPSFromDistMatrixFile(ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&,bool);
		// rips complex from delay reconstruction of time series
		bool buildBRIPSFromTimeSeriesFile(ifstream&, MComplex<CC,BC>&, RIPS<CC,double,BC>&);
		
		//	Compute barcode
		void ComputeBarcode(string input_type, string infile, string outfile, string eng);
	
};
