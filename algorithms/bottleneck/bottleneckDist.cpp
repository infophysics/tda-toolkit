//============================================================================
// Name        : bottleneckDist.cpp
// Author      : Jesse Berwald (modified Miro's code)
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <iostream>
#include <fstream>
#include "CBottleneckDistance.cpp"

using namespace std;


int main(int argc, char* argv[])
{
    // must pass in two args, one for each persdia file, so #argc = 1 + 2
    // may also want to pass in max_gen later
    if (argc < 3 || argc > 4) {
        cout
            << endl
            << "Must pass in filename1 and filename2"
            << endl
            << "to find the distance between diagrams."
            << endl;
        return -1;
    }

    // file names
    char* file1 = argv[1];
    char* file2 = argv[2];

    CBottleneckDistance d;

    double max_gen = 1;

    // compute the distance between persistence diagrams
    double distance = d.Distance(file1, file2, max_gen);

    cout << distance << endl;

    return 0;
}
