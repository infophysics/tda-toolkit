/**
 * BottleneckDistanceDouble.cpp
 */

#include <iostream>
#include <fstream>

#include "CBottleneckDistance.h"
using namespace std;


int main() {

    CBottleneckDistance d;

    int n_files = 7;
    int B = 1;

    double max_gen = 1;

    //  double skip_last = 1;
    //  double skip_short = 0;
    //  d.SkipLastGenerators( 4.97 );
    //  d.SkipShortGenerators( 0.05 );

    double *distances;

    distances = new double[ n_files * n_files];

    for (int i = 0; i < n_files; ++i)
        distances[ i * n_files ] = 0;

    for (int file_1 = 0; file_1 < n_files; ++file_1)
        for (int file_2 = file_1 + 1; file_2 < n_files; ++file_2) {
            char file_1_name[100];
            char file_2_name[100];
            sprintf (file_1_name, "./Inv_file_%d_B%d.txt", file_1, B);
            sprintf (file_2_name, "./Inv_file_%d_B%d.txt", file_2, B);

            double distance = d.Distance( file_1_name, file_2_name, max_gen);
            distances[ file_1 * n_files + file_2] = distance;
            distances[ file_2 * n_files + file_1] = distance;

            std::cout << "Distance[ " << file_1 << " , " << file_2 << " ] = " << distance << "\n";
        }

    char output_file_name[100];
    sprintf ( output_file_name, "./Inv_Bottleneck_B%d.txt", B);
    std::ofstream out_file( output_file_name );
    if (!out_file){
        std::cout<< "Unable to open input file " <<  output_file_name << "\n";
        return false;
    }

    for (int i = 0; i < n_files; ++i) {
        for (int j = 0; j < n_files; ++j) {
            out_file << distances[ i * n_files + j ] << " ";
        }
        out_file << "0 \n";
    }

    for (int i = 0; i < n_files; ++i)
        out_file << "0 ";

    out_file << "0 \n";
    out_file.close();

    delete distances;

    return 0;
}

