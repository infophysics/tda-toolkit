/**
 * Author: Kelvin Abrokwa (kelvinabrokwa@gmail.com)
 */

#include "mex.h"
#include "CBottleneckDistance.cpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void printGenerators(std::vector<Generator> generators);
std::vector<Generator> loadGenerators(double* A, int m, int n, double maxLevel);


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *A, *B;
    int mA, nA, mB, nB;

    // TODO: accept as parameter
    double maxLevel = 1;

    // validate Matlab arguments
    if (nrhs != 2)
        mexErrMsgIdAndTxt("Bottleneck:b:nrhs", "Arguments: 2 Mx2 matrices");

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    mA = mxGetM(prhs[0]);
    nA = mxGetN(prhs[0]);
    mB = mxGetM(prhs[1]);
    nB = mxGetN(prhs[1]);

    if (nA != 2 || nB != 2)
        mexErrMsgIdAndTxt("Bottleneck:b:nrhs", "Mx2 matrices required");

    CBottleneckDistance d;

    std::vector<Generator> gensA = loadGenerators(A, mA, nB, maxLevel);
    std::vector<Generator> gensB = loadGenerators(B, mB, nB, maxLevel);

    float distance = d.Distance(gensA, gensB, maxLevel);

    plhs[0] = mxCreateDoubleScalar(distance);
}


/**
 *
 */
std::vector<Generator> loadGenerators(double* A, int m, int n, double maxLevel)
{
    std::vector<Generator> generators;

    Generator gen;

    for (int row = 0; row < m; row++) {
        for (int col = 0; col < n; col++) {
            if (col % 2)
                gen.death = A[row + col*m] == -1 ? maxLevel : A[row + col*m];
            else
                gen.birth = A[row + col*m] == -1 ? maxLevel : A[row + col*m];
        }
        // TODO: incorporate generator neglection
        generators.push_back(gen);
    }

    return generators;
}




