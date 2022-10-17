/*=========================================================================
 *
 * 09-12-2020 - coded by Alessandro Muscoloni
 *
 * > INPUT
 * x - adjacency matrix (sparse, double, symmetric, zero-diagonal)
 * G - matrix of geodesics (full, double, symmetric, zero-diagonal)
 * T - matrix of topological shortest paths (full, double, symmetric, zero-diagonal)
 * order - vector indicating the order in which the nodes should be processed
 *
 * > OUTPUT
 * P - matrix indicating for each node pair (i,j) the mean geometrical projection of topological shortest paths from i to j
 *
 * Usage:
 * P = compute_pTSP_mex(x, G, T, order)
 *
 * In Windows, install the MATLAB Support for MinGW-w64 C/C++ Compiler and then compile with:
 * mex C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a compute_pTSP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *
 *=======================================================================*/

#include "mex.h"
#include "matrix.h"
#include <omp.h>

void compute_tsp_proj_rec(mwSize N, mwIndex *ir, mwIndex *jc, mxDouble *G, mxDouble *T, mxDouble *P,
                          mxDouble *L, mwSize *Pcount, mxLogical *inpath, mwIndex s, mwSize M, mwIndex m,
                          mwIndex u, mwSize lt, mxDouble lg, mxDouble g)
{
    /* add the node to the path */
    inpath[u*M+m] = true;
    lt = lt + 1;
    lg = lg + g;
    
    /* store the projection when needed */
    if ((lt==T[u*N+s]) && ((L[s]>L[u]) || ((L[s]==L[u]) && (s<u))))
    {
        Pcount[u*M+m] = Pcount[u*M+m] + 1;
        P[u*N+s] = P[u*N+s] + lg;
    }
    
    /* if maximum path length not reached, continue recursion
     * visiting the neighbours not already in the path */
    if (lt<L[s])
    {
        for (mwIndex i=jc[u]; i<jc[u+1]; i++)
        {
            if (!inpath[ir[i]*M+m])
            {
                compute_tsp_proj_rec(N, ir, jc, G, T, P, L, Pcount, inpath, s, M, m, ir[i], lt, lg, G[ir[i]*N+u]);
            }
        }
    }

    /* remove the node from the path */
    inpath[u*M+m] = false;
}

/*=======================================================================*/

void compute_tsp_proj_main(mwSize N, mwIndex *ir, mwIndex *jc, mxDouble *G, mxDouble *T, mxDouble *P,
                           mxDouble *L, mwSize *Pcount, mxLogical *inpath, mwIndex s, mwSize M, mwIndex m)
{
    /* initialize variables */
    mwSize lt = 0;
    mxDouble lg = 0;
    mwIndex i;
    for (i=0; i<N; i++)
    {
        Pcount[i*M+m] = 0;
        inpath[i*M+m] = false;
    }
    inpath[s*M+m] = true;
    
    /* for each neighbour, start recursive computation */
    for (i=jc[s]; i<jc[s+1]; i++)
    {
        compute_tsp_proj_rec(N, ir, jc, G, T, P, L, Pcount, inpath, s, M, m, ir[i], lt, lg, G[ir[i]*N+s]);
    }
    
    /* update output */
    for (i=0; i<N; i++)
    {
        if ((L[s]>L[i]) || ((L[s]==L[i]) && (s<i)))
        {
            P[i*N+s] = P[i*N+s] / (mxDouble) Pcount[i*M+m];
            P[s*N+i] = P[i*N+s];
        }
    }
}

/*=======================================================================*/

void compute_L(mwSize N, mxDouble *T, mxDouble *order, mxDouble *L)
{
    mxLogical *mask = (mxLogical*) mxCalloc(N, sizeof(mxLogical));
    mwIndex i, j, s;
    for (i=0; i<N; i++)
    {
        s = (mwIndex) order[i]-1;
        for (j=0; j<N; j++)
        {
            if ((mask[j]==false) && (T[j*N+s] > L[s]))
            {
                L[s] = T[j*N+s];
            }
        }
        mask[s] = true;
    }
    mxFree(mask);
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
    mxDouble *T = mxGetPr(prhs[2]);
    mxDouble *order = mxGetPr(prhs[3]);
    
    /* allocate output */
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    mxDouble *P = mxGetPr(plhs[0]);
    
    /* compute L: for each node indicates the maximum path length to compute */
    mxDouble *L = (mxDouble*) mxCalloc(N, sizeof(mxDouble));
    compute_L(N, T, order, L);
    
    /* allocate shared support variables */
    mwSize M = omp_get_max_threads();
    mwSize *Pcount = (mwSize*) mxCalloc(M*N, sizeof(mwSize));
    mxLogical *inpath = (mxLogical*) mxCalloc(M*N, sizeof(mxLogical));
    
    /* for each source node, run computational routine */
    #pragma omp parallel for shared(N, ir, jc, G, T, order, P, L, Pcount, inpath, M) schedule(dynamic)
    for (mwIndex i=0; i<N-1; i++)
    {
        mwIndex s = (mwIndex) order[i]-1;
        mwIndex m = (mwIndex) omp_get_thread_num();
        compute_tsp_proj_main(N, ir, jc, G, T, P, L, Pcount, inpath, s, M, m);
    }
    
    /* free memory */
    mxFree(L);
    mxFree(Pcount);
    mxFree(inpath);
}
