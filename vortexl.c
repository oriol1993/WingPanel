/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void arrayProduct(double *p, double *v1, double *v2, double *v)
{
    double n1[3], n2[3], cp[3], r1, r2, cp2, vn[3], r0r1, r0r2, K;
    
    n1[0] = p[0]-v1[0];
    n1[1] = p[1]-v1[1];
    n1[2] = p[2]-v1[2];
    n2[0] = p[0]-v2[0];
    n2[1] = p[1]-v2[1];
    n2[2] = p[2]-v2[2];    
    cp[0] = n1[1]*n2[2]-n1[2]*n2[1];
    cp[1] = n1[2]*n2[0]-n1[0]*n2[2];
    cp[2] = n1[0]*n2[1]-n1[1]*n2[0];
    r1 = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
    r2 = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
    cp2 = cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2];
    vn[0] = v2[0]-v1[0];
    vn[1] = v2[1]-v1[1];
    vn[2] = v2[2]-v1[2];
    r0r1 = vn[0]*n1[0]+vn[1]*n1[1]+vn[2]*n1[2];
    r0r2 = vn[0]*n2[0]+vn[1]*n2[1]+vn[2]*n2[2];
    K = (r0r1/r1-r0r2/r2)/(4*3.14159265359*cp2);
    v[0] = cp[0]*K;
    v[1] = cp[1]*K;
    v[2] = cp[2]*K;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *p, *v1, *v2;               /* 1xN input matrix */
    double *v;              /* output matrix */

    /* create a pointer to the real data in the input matrix  */
    p = mxGetPr(prhs[0]);
    v1 = mxGetPr(prhs[1]);
    v2 = mxGetPr(prhs[2]);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,3,mxREAL);
    /* get a pointer to the real data in the output matrix */
    v = mxGetPr(plhs[0]);

    /* call the computational routine */
    arrayProduct(p,v1,v2,v);
}
