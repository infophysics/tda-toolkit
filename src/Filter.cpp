/*
 * 	PH filter functions for various data types				Nick Carrara
 */
#include "Filter.h"
using namespace std;

//	Filter 2D
Filter2D::Filter2D(){}
Filter2D::~Filter2D(){}

void Filter2D::loadBinaryFromFile(const char* input_file, std::string format){
	
	if(format == "DIPHA" || format == "dipha"){ 
		ifstream reading_file; 
		if (input_file && file_stream.fail()) {
			cerr << "ERROR! Could not open file " << filename << std::endl;
			return;
		}
		ifstream fin( filename, ios::in | ios::binary ); 
		int64_t d; double ax; double ay;
		fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
		assert(d == 8067171840);
		fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
		assert(d == 1);
		fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
		fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim
		assert(d == 2);
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		ax = d;
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		ay = d;
		assert(0 < ax && 0 < ay);
	
		double dense[ax][ay];
		double dou;
		for (int y = 0; y < ay + 2; ++y) {
			for (int x = 0; x < ax + 2; ++x) {
				if(0 < x && x <= ax && 0 < y && y <= ay){
					if (!fin.eof()) {
						fin.read( ( char * ) &dou, sizeof( double ) );
						dense[x][y] = dou;
					} else {
						cout << "ERROR! End of file!" << endl;
					}
				}
				else {
					dense[x][y] = threshold;
				}
			}
		}
		fin.close();
		for (int x = 0; x < ax; x++){
			vector<int> points;
			for (int y = 0; y < ay; y++){
				points.push_back(dense[x][y]);
			}
			m_Binary.push_back(points);
		}
	}  else if(format == "PERSEUS" || format == "perseus"){
		ifstream reading_file; 
		reading_file.open(filename.c_str(), ios::in);
		int64_t d; double ax; double ay;
		string reading_line_buffer; 
		getline(reading_file, reading_line_buffer); 
		d = atoi(reading_line_buffer.c_str()); 
		getline(reading_file, reading_line_buffer); 
		ax = atoi(reading_line_buffer.c_str()); 
		getline(reading_file, reading_line_buffer); 
		ay = atoi(reading_line_buffer.c_str()); 
		assert(0 < ax && 0 < ay);
		double dense[ax][ay];
		for (int y = 0; y <ay + 2; ++y) { 
			for (int x = 0; x < ax + 2; ++x) { 
				if(0 < x && x <= ax && 0 < y && y <= ay){ 
					if (!reading_file.eof()) { 
						getline(reading_file, reading_line_buffer); 
						dense[x][y] = atoi(reading_line_buffer.c_str()); 
						if (dense[x][y] == -1) { 
							dense[x][y] = 99999;
						} 
					} 
				} 
				else { 
					dense[x][y] = 99999; 
				} 
			} 
		}
		for (int x = 0; x < ax; x++){
			vector<int> points;
			for (int y = 0; y < ay; y++){
				points.push_back(dense[x][y]);
			}
			m_Binary.push_back(points);
		}
	} 
	else{
		ifstream file(input_file);
	 
		vector<vector<string> > feed;
	 
		string line = "";
		while (getline(file, line))
		{
			vector<string> vec;
			boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
			feed.push_back(vec);
		}
		for (int x = 0; x < feed.size(); x++){
			vector<int> points;
			for (int y = 0; y < feed[x].size(); y++){
				points.push_back(feed[x][y]);
			}
			m_Binary.push_back(points);
		}
	}
}
void Filter2D::loadContinuousFromFile(const char* input_file, std::string format){

}

void Filter2D::countDeadCells(){
	int deadCells;
	for (int i = 0; i < m_Binary.size(); i++){
		for (int j = 0; j < m_Binary[i].size(); j++){
			if (m_Binary[i][j] == 0) deadCells++;
		}
	}
	m_DeadCells = deadCells;
}
		
//	Various filterings
//	Binary filterings
void Filter2D::filterBinaryVonNeumann(int threshold){	//	only hard boundary right now
	//	count the number of dead cells
	int deadCells = m_DeadCells;
	bool thres = true;
	int currentState = 1;
	//	loop
	while(deadCells > 0 || thres == false){
		vector<vector<int> > binaryCopy = m_Binary;
		for (int i = 0; i < binaryCopy.size(); i++){
			for (int j = 0; j < binaryCopy[i].size(); j++){
				//	Check state
				if (binaryCopy[i][j] == currentState){
					//	Look at all neighbors
					if (i == 0){
						if (j == 0){
							if (binaryCopy[0][1] == 0) binaryCopy[0][1] = currentState + 1; deadCells--;
							if (binaryCopy[1][0] == 0) binaryCopy[1][0] = currentState + 1; deadCells--;
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[0][j-2] == 0) binaryCopy[0][j-2] = currentState + 1;deadCells--;
							if (binaryCopy[1][j-1] == 0) binaryCopy[1][j-1] = currentState + 1;deadCells--;
						}
						else{
							if (binaryCopy[0][j-1] == 0) binaryCopy[0][j-1] = currentState + 1;deadCells--;
							if (binaryCopy[0][j+1] == 0) binaryCopy[1][j+1] = currentState + 1;	deadCells--;
							if (binaryCopy[1][j] == 0) binaryCopy[1][j] = currentState + 1;deadCells--;
						}
					}
					else if (i == binaryCopy[i].size() - 1){
						if (j == 0){
							if (binaryCopy[i][1] == 0) binaryCopy[i][1] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][0] == 0) binaryCopy[i-1][0] = currentState + 1;deadCells--;
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-2] == 0) binaryCopy[i][j-2] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][j-1] == 0) binaryCopy[i-1][j-1] = currentState + 1;deadCells--;
						}
						else{
							if (binaryCopy[i][j-1] == 0) binaryCopy[i][j-1] = currentState + 1;deadCells--;
							if (binaryCopy[i][j+1] == 0) binaryCopy[i][j+1] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][j] == 0) binaryCopy[i-1][j] = currentState + 1;deadCells--;
						}	
					}
					else{
						if (j == 0){
							if (binaryCopy[i][1] == 0) binaryCopy[i][1] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][0] == 0) binaryCopy[i-1][0] = currentState + 1;deadCells--;
							if (binaryCopy[i+1][0] == 0) binaryCopy[i+1][0] = currentState + 1;deadCells--;
						}
						else if (j == binaryCopy[i].size()-1){
							if (binaryCopy[i][j-2] == 0) binaryCopy[i][j-2] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][j-1] == 0) binaryCopy[i-1][j-1] = currentState + 1;deadCells--;
							if (binaryCopy[i+1][j-1] == 0) binaryCopy[i+1][j-1] = currentState + 1;deadCells--;
						}
						else{
							if (binaryCopy[i][j-1] == 0) binaryCopy[i][j-1] = currentState + 1;deadCells--;
							if (binaryCopy[i][j+1] == 0) binaryCopy[i][j+1] = currentState + 1;deadCells--;
							if (binaryCopy[i-1][j] == 0) binaryCopy[i-1][j] = currentState + 1;deadCells--;
							if (binaryCopy[i+1][j] == 0) binaryCopy[i+1][j] = currentState + 1;deadCells--;
						}	
					}
				}
			}
		}
		if (currentState >= threshold) thres = false;
	}
	m_Binary = binaryCopy;
}
void Filter2D::filterBinaryMoore(int threshold){
	
}


