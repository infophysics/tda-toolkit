/*

 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License
 
 Copyright (c) 2015–2018 Ulrich Bauer
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/


//#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

#include <string>

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;

enum file_format3 {
	LOWER_DISTANCE_MATRIX,
	UPPER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	DIPHA3,
	RIPSER
};

class Ripser{
	private:
		std::vector<std::vector<float> > m_barcode;
	public:
		Ripser();
		virtual ~Ripser();
		void print_usage_and_exit(int exit_code);
		void ComputeBarcode(const char* filename, long dim, float thres, float rat, std::string form, long mod);
		
		std::vector<std::vector<float> > getBarcode(){ return m_barcode; }
		void saveBarcodeToFile(std::string output_file);
};







