#pragma once
/*
 * 	PH filter functions for various data types				Nick Carrara
 */
#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <sstream>
#include <cstdint>
#include <math.h>

class Filter2D{
	
	private:
		std::vector<std::vector<double> > m_Binary;
		std::vector<std::vector<double> > m_Continuous;
		int m_DeadCells;
		
	public:
		Filter2D();
		virtual ~Filter2D();
		void loadBinaryFromFile(const char* input_file);
		void loadContinuousFromFile(const char* input_file);
		void countDeadCells();
		void printFiltration(std::vector<std::vector<int> > &binaryCopy);
		//	Various filterings
		//	Binary filterings
		void filterBinaryL1(double threshold);
		void filterBinaryL2(double threshold);
		void filterBinaryLinf(double threshold);
		void filter3StateAsBinary(double alive);		
		//	Save filtration
		void saveBinaryFiltration(const char* output_file);
};
