/*=========================================================================
 *
 * 22-12-2020 - coded by Alessandro Muscoloni
 *
 * > INPUT
 * x - adjacency matrix (sparse, double, symmetric, zero-diagonal)
 * G - matrix of geodesics (full, double, symmetric, zero-diagonal)
 *
 * > OUTPUT
 * P - matrix indicating for each node pair (i,j) the geometrical path length of GR from i to j
 *
 * Usage:
 * P = compute_pGRP_mex(x, G)
 *
 * In Windows, install the MATLAB Support for MinGW-w64 C/C++ Compiler and then compile with:
 * mex C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a compute_pGRP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *
 *=======================================================================*/

#include "mex.h"
#include "matrix.h"
#include <omp.h>

void greedy_routing_to_j_rec(mwSize N, mwIndex *ir, mwIndex *jc, mxDouble *G, mxDouble *P,
                             mwIndex *next, mwIndex m, mwIndex j, mwIndex i)
{
    mwIndex k = next[m*N+i];
    if (k==j)
    {
        P[j*N+i] = G[j*N+i];
    }
    else if (next[m*N+k]==i)
    {
        P[j*N+i] = mxGetInf();
    }
    else
    {
        if (P[j*N+k]==0)
        {
            greedy_routing_to_j_rec(N, ir, jc, G, P, next, m, j, k);
        }
        if (P[j*N+k]==mxGetInf())
        {
            P[j*N+i] = mxGetInf();
        }
        else
        {
            P[j*N+i] = G[k*N+i] + P[j*N+k];
        }
    }   
}

/*=======================================================================*/

void greedy_routing_to_j(mwSize N, mwIndex *ir, mwIndex *jc, mxDouble *G, mxDouble *P,
                         mwIndex *next, mwIndex m, mwIndex j)
{
    mwIndex i;
    for (i=0; i<N; i++)
    {
        if ((i!=j) && (P[j*N+i]==0))
        {
            greedy_routing_to_j_rec(N, ir, jc, G, P, next, m, j, i);
        }
    }
}

/*=======================================================================*/

void find_next_to_j(mwSize N, mwIndex *ir, mwIndex *jc, mxDouble *G,
                    mwIndex *next, mwIndex m, mwIndex j)
{
    mwIndex i, k;
    mxDouble Gmin;
    for (i=0; i<N; i++)
    {
        next[m*N+i] = -1;
        if (i==j)
        {
            continue;
        }
        Gmin = mxGetInf();
        for (k=jc[i]; k<jc[i+1]; k++)
        {
            if (ir[k]==j)
            {
                next[m*N+i] = j;
                break;
            }
            else if (G[j*N+ir[k]] < Gmin)
            {
                next[m*N+i] = ir[k];
                Gmin = G[j*N+ir[k]];
            }
        }
    }
}

/*=======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* get input */
    mwSize N = mxGetN(prhs[0]);
    mwIndex *ir = mxGetIr(prhs[0]);
    mwIndex *jc = mxGetJc(prhs[0]);
    mxDouble *G = mxGetPr(prhs[1]);
    
    /* allocate output */
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    mxDouble *P = mxGetPr(plhs[0]);
    
    /* allocate shared support variable */
    mwSize M = omp_get_max_threads();
    mwIndex *next = (mwIndex*) mxCalloc(N*M, sizeof(mwIndex));
    
    /* for each destination node, run computational routine */
    #pragma omp parallel for shared(N, ir, jc, G, P, next) schedule(dynamic)
    for (mwIndex j=0; j<N; j++)
    {
        mwIndex m = (mwIndex) omp_get_thread_num();
        find_next_to_j(N, ir, jc, G, next, m, j);
        greedy_routing_to_j(N, ir, jc, G, P, next, m, j);
    }
    
    /* free memory */
    mxFree(next);
}
