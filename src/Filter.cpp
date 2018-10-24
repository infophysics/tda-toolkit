/*
 * 	PH filter functions for various data types				Nick Carrara
 */
#include "Filter.h"
using namespace std;

//	Filter 2D
Filter2D::Filter2D(){}
Filter2D::~Filter2D(){}

void Filter2D::loadBinaryFromFile(const char* input_file){
		std::ifstream inputFile(input_file);
		int l = 0;
		 
		while (inputFile) {
			l++;
			std::string s;
			if (!std::getline(inputFile, s)) break;
		        if (s[0] != '#') {
		            std::istringstream ss(s);
		            std::vector<double> record;
		 
		            while (ss) {
		                std::string line;
		                if (!std::getline(ss, line, ','))
		                    break;
		                try {
		                    record.push_back(stof(line));
		                }
		                catch (const std::invalid_argument e) {
		                    std::cout << "NaN found in file " << input_file << " line " << l
		                         << std::endl;
		                    e.what();
		                }
		            }
		 
		            m_Binary.push_back(record);
		       }
		   }
		 
		   if (!inputFile.eof()) {
		       std::cerr << "Could not read file " << input_file << "\n";
		       std::__throw_invalid_argument("File not found.");
		   }
		
}
void Filter2D::loadContinuousFromFile(const char* input_file){

}

void Filter2D::countDeadCells(){
	int deadCells = 0;
	for (int i = 0; i < m_Binary.size(); i++){
		for (int j = 0; j < m_Binary[i].size(); j++){
			if (m_Binary[i][j] == 0) deadCells++;
		}
	}
	m_DeadCells = deadCells;
}
		
//	Various filterings
//	Binary filterings
void Filter2D::filterBinaryL1(double threshold){	//	only hard boundary right now
	//	count the number of dead cells
	countDeadCells();
	int deadCells = m_DeadCells;
	bool thres = false;
	double currentState = 1.0;
	//	loop
	vector<vector<double> > binaryCopy = m_Binary;
	while(deadCells > 0 && thres == false){
		for (int i = 0; i < binaryCopy.size(); i++){
			for (int j = 0; j < binaryCopy[i].size(); j++){
				//	Check state
				if (binaryCopy[i][j] == currentState){
					//	Look at all neighbors
					if (i == 0){
						if (j == 0){
							if (binaryCopy[0][1] == 0){
								binaryCopy[0][1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[1][0] == 0){
								binaryCopy[1][0] = currentState + 1; 
								deadCells--;
							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[0][j-1] == 0){
								binaryCopy[0][j-1] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[1][j] == 0){
								binaryCopy[1][j] = currentState + 1;
								deadCells--;
							}
						}
						else{
							if (binaryCopy[0][j-1] == 0){
								binaryCopy[0][j-1] = currentState + 1;
								deadCells--;		
							}
							if (binaryCopy[0][j+1] == 0){
								binaryCopy[0][j+1] = currentState + 1;	
								deadCells--;
							}
							if (binaryCopy[1][j] == 0){
								binaryCopy[1][j] = currentState + 1;
								deadCells--;
							}
						}
					}
					else if (i == binaryCopy.size() - 1){
						if (j == 0){
							if (binaryCopy[i][1] == 0){
								binaryCopy[i][1] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[i-1][0] == 0){
								binaryCopy[i-1][0] = currentState + 1;
								deadCells--;
							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;		

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
						}
						else{
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i][j+1] == 0){
								binaryCopy[i][j+1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
						}	
					}
					else{
						if (j == 0){
							if (binaryCopy[i][1] == 0){
								binaryCopy[i][1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][0] == 0){
								binaryCopy[i-1][0] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][0] == 0){
								binaryCopy[i+1][0] = currentState + 1;
								deadCells--;

							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][j] == 0){
								binaryCopy[i+1][j] = currentState + 1;
								deadCells--;

							}
						}
						else{
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i][j+1] == 0){
								binaryCopy[i][j+1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][j] == 0){
								binaryCopy[i+1][j] = currentState + 1;
								deadCells--;
							}
						}	
					}
				}
			}
		}
		currentState++;
		if (currentState >= threshold) thres = true;
	}
	m_Binary = binaryCopy;
}

void Filter2D::filterBinaryL2(double threshold){
	//	Make a list of all the 1 cells
	std::vector<std::vector<double> > liveCells;
	for (int i = 0; i < m_Binary.size(); i++){
		for (int j = 0; j < m_Binary[i].size(); j++){
			if (m_Binary[i][j] == 1) liveCells.push_back({i,j});
		}
	}
	//	Now determine the radius for each zero
	for (int i = 0; i < m_Binary.size(); i++){
		for (int j = 0; j < m_Binary[i].size(); j++){
			if (m_Binary[i][j] == 0){
				double min_distance = 10.0e10;
				for (int k = 0; k < liveCells.size(); k++){
					double temp_distance = (i - liveCells[k][0])*(i - liveCells[k][0]) + (j - liveCells[k][1])*(j - liveCells[k][1]);
					temp_distance = sqrt(temp_distance);
					if (temp_distance <= min_distance) min_distance = temp_distance;
				}
				m_Binary[i][j] = min_distance + 1;
			}
		}
	}
}




void Filter2D::filterBinaryLinf(double threshold){	//	only hard boundary right now
	//	count the number of dead cells
	countDeadCells();
	int deadCells = m_DeadCells;
	bool thres = false;
	int currentState = 1;
	//	loop
	vector<vector<double> > binaryCopy = m_Binary;
	while(deadCells > 0 && thres == false){
		for (int i = 0; i < binaryCopy.size(); i++){
			for (int j = 0; j < binaryCopy[i].size(); j++){
				//	Check state
				if (binaryCopy[i][j] == currentState){
					//	Look at all neighbors
					if (i == 0){
						if (j == 0){
							if (binaryCopy[0][1] == 0){
								binaryCopy[0][1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[1][0] == 0){
								binaryCopy[1][0] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[1][1] == 0){
								binaryCopy[1][1] = currentState + 1; 
								deadCells--;
							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[0][j-1] == 0){
								binaryCopy[0][j-1] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[1][j] == 0){
								binaryCopy[1][j] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[1][j-1] == 0){
								binaryCopy[1][j-1] = currentState + 1; 
								deadCells--;
							}
						}
						else{
							if (binaryCopy[0][j-1] == 0){
								binaryCopy[0][j-1] = currentState + 1;
								deadCells--;		
							}
							if (binaryCopy[0][j+1] == 0){
								binaryCopy[0][j+1] = currentState + 1;	
								deadCells--;
							}
							if (binaryCopy[1][j] == 0){
								binaryCopy[1][j] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[1][j-1] == 0){
								binaryCopy[1][j-1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[1][j+1] == 0){
								binaryCopy[1][j+1] = currentState + 1; 
								deadCells--;
							}
						}
					}
					else if (i == binaryCopy.size() - 1){
						if (j == 0){
							if (binaryCopy[i][1] == 0){
								binaryCopy[i][1] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[i-1][0] == 0){
								binaryCopy[i-1][0] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[i-1][1] == 0){
								binaryCopy[i-1][1] = currentState + 1; 
								deadCells--;
							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;		

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j-1] == 0){
								binaryCopy[i-1][j-1] = currentState + 1; 
								deadCells--;
							}
						}
						else{
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i][j+1] == 0){
								binaryCopy[i][j+1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j-1] == 0){
								binaryCopy[i-1][j-1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i-1][j+1] == 0){
								binaryCopy[i-1][j+1] = currentState + 1; 
								deadCells--;
							}
						}	
					}
					else{
						if (j == 0){
							if (binaryCopy[i][1] == 0){
								binaryCopy[i][1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][0] == 0){
								binaryCopy[i-1][0] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][0] == 0){
								binaryCopy[i+1][0] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][1] == 0){
								binaryCopy[i-1][1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i+1][1] == 0){
								binaryCopy[i+1][1] = currentState + 1; 
								deadCells--;
							}
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][j] == 0){
								binaryCopy[i+1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j-1] == 0){
								binaryCopy[i-1][j-1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i+1][j-1] == 0){
								binaryCopy[i+1][j-1] = currentState + 1; 
								deadCells--;
							}
						}
						else{
							if (binaryCopy[i][j-1] == 0){
								binaryCopy[i][j-1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i][j+1] == 0){
								binaryCopy[i][j+1] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i-1][j] == 0){
								binaryCopy[i-1][j] = currentState + 1;
								deadCells--;

							}
							if (binaryCopy[i+1][j] == 0){
								binaryCopy[i+1][j] = currentState + 1;
								deadCells--;
							}
							if (binaryCopy[i-1][j-1] == 0){
								binaryCopy[i-1][j-1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i-1][j+1] == 0){
								binaryCopy[i-1][j+1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i+1][j-1] == 0){
								binaryCopy[i+1][j-1] = currentState + 1; 
								deadCells--;
							}
							if (binaryCopy[i+1][j+1] == 0){
								binaryCopy[i+1][j+1] = currentState + 1; 
								deadCells--;
							}
						}	
					}
				}
			}
		}
		currentState++;
		if (currentState >= threshold) thres = true;
	}
	m_Binary = binaryCopy;
}

void Filter2D::filter3StateAsBinary(double alive){
	for (int i = 0; i < m_Binary.size(); i++){
		for (int j = 0; j < m_Binary[i].size(); j++){
			if (m_Binary[i][j] == alive) m_Binary[i][j] = 1;
			else m_Binary[i][j] = 0;
		}
	}
}


void Filter2D::printFiltration(std::vector<std::vector<int> > &binaryCopy){
	for (int i = 0; i < binaryCopy.size(); i++){
		for (int j = 0; j < binaryCopy[i].size(); j++){
			std::cout << binaryCopy[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}


void Filter2D::saveBinaryFiltration(const char* output_file){
	std::ofstream myfile;
		myfile.open(output_file);
		for (int i = 0; i < m_Binary.size(); i++){
			for (int j = 0; j < m_Binary[i].size(); j++){
				myfile << m_Binary[i][j];
				if (j < m_Binary[i].size() - 1){
					myfile << ",";
				}
			}
			myfile << "\n";
		}
		myfile.close();
	
}

