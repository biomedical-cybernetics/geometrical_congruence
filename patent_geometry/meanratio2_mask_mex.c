/*=========================================================================
 *
 * 22-12-2020 - coded by Alessandro Muscoloni
 *
 * > INPUT
 * num1 - matrix indicating the first numerator (full, double)
 * num2 - matrix indicating the second numerator (full, double)
 * den - matrix indicating the denominator (full, double)
 * mask - matrix indicating the elements to consider (full, logical)
 *
 * > OUTPUT
 * mr1 - equivalent to mean(num1(mask)./den(mask))
 * mr2 - equivalent to mean(num2(mask)./den(mask))
 *
 * Usage:
 * [mr1, mr2] = meanratio2_mask_mex(num1, num2, den, mask)
 *
 * In Windows, install the MATLAB Support for MinGW-w64 C/C++ Compiler and then compile with:
 * mex meanratio2_mask_mex.c
 *
 *=======================================================================*/

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* get input */
    mwSize N = mxGetN(prhs[0]);
    mxDouble *num1 = mxGetPr(prhs[0]);
    mxDouble *num2 = mxGetPr(prhs[1]);
    mxDouble *den = mxGetPr(prhs[2]);
    mxLogical *mask = mxGetLogicals(prhs[3]);
    
    /* run computational routine */
    mxDouble c = 0;
    mxDouble s1 = 0;
    mxDouble s2 = 0;
    mwIndex i, j;
    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            if (mask[j*N+i])
            {
                c = c + 1;
                s1 = s1 + (num1[j*N+i]/den[j*N+i]);
                s2 = s2 + (num2[j*N+i]/den[j*N+i]);
            }
        }
    }
        
    /* allocate output */
    plhs[0] = mxCreateDoubleScalar(s1/c);
    plhs[1] = mxCreateDoubleScalar(s2/c);
}
