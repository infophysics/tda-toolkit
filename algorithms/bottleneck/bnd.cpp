/**
 * Author: Kelvin Abrokwa (kelvinabrokwa@gmail.com)
 *
 * Bottle Neck Distance Wrapper for Matlab mex
 * Takes 2 string arguments in Matlab
 */
#include "mex.h"
#include "CBottleneckDistance.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    // validate Matlab input
    if (nrhs != 2)
        mexErrMsgIdAndTxt("MyToolbox:bn:nrhs", "Need 2 arguments");

    // read arguments
    char *file1 = mxArrayToString(prhs[0]);
    char *file2 = mxArrayToString(prhs[1]);

    double maxGen = 1;

    CBottleneckDistance d;

    // compute the distance between two persistance diagrams
    double distance = d.Distance( file1, file2, maxGen);

    // return the distance to the Matlab caller
    plhs[0] = mxCreateDoubleScalar(distance);
}
