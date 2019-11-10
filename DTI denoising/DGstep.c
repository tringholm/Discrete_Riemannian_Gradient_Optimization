// dgstepMex propagates one timestep of the discrete gradient algorithm 
// using a C implementation
// 
// Input:
// u       - current restored image
// Nx & Ny - image dimensions
// g       - noisy input greyscale image
// dt      - time step size
// a       - total variation regularization weight
// b       - curvature term regularization weight
// s       - fidelity term is computed in L^s norm
// epsilon - smoothing parameter
// tol     - tolerance for Brent-Dekker nonlinear equation solver
// K       - logical map with false values on pixels to be inpainted
// 
// Output:
// u       - updated image after one time step
// fcount  - total number of function used in timestep
// 
// Torbjørn Ringholm
// Email           : torbjorn.ringholm@ntnu.no
// Last updated    : 20/11/2017

#define PI 3.14159265358979323846

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <matrix.h>
#include "mex.h"

#define g_in            prhs[0]
#define u_in            prhs[1]
#define mask_in         prhs[2]
#define dp_in           prhs[3]
#define dqx_in          prhs[4]
#define dqy_in          prhs[5]
#define Suinvprod_in    prhs[6]
#define Sginvprod_in    prhs[7]
#define p_in            prhs[8]
#define gamma_in        prhs[9]
#define tol_in          prhs[10]
#define xtol_in         prhs[11]
#define dt_in           prhs[12]
#define M_in            prhs[13]
#define N_in            prhs[14]
#define Tmax_in         prhs[15]

// ---------------------------------------------- mex version of eigens for 3x3 matrix
void eignvec3x3(double *V, double *L, double *out);

// ---------------------------------------------- mex version of product simplifier
void prodform(double *V, double *L, double *out);


// ---------------------------------------------- mex version of inverse 3x3 matrix
void inv3x3(double *uT, double* uTinv);

// ---------------------------------------------- mex version of main loop
void DGstep(double *g, double *u, int *mask, double *dp, double *dqx,
        double *dqy, double *Suinvprod, double *Sginvprod, int p, double gamma,
        double tol, double xtol, double dt, int M, int N, int Tmax, double *residual, double *Vhist);

// ---------------------------------------------- mex version of Brent-Dekker
void BrentDekker(double *x, double *uT, double *uTinv, double *uxp, double *uyp,
        double *uxm, double *uym, double *g, int p, double gamma, int l, double Vo, double dt, int *masks,
        double *out, double* dug, double *dxp, double *dxm, double* dyp, double* dym, double*un, double tol, int interval);

// ---------------------------------------------- mex version of coordinate evaluation
void Vdiff_solv(double alpha, double *uT, double *uTinv, double *uxp, double *uyp,
        double *uxm, double *uym, double *g, int p, double gamma, int l, double Vo, double dt, int *masks,
        double *out, double* dug, double *dxp, double *dxm, double* dyp, double* dym, double* un);

// ---------------------------------------------- mex version of manifold distance
double d_M(double *XX, double *A);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------ initialize and get input
    double gamma, tol, xtol, dt;
    double *g, *u, *dp, *dqx, *dqy, *Suinvprod, *Sginvprod, *residual, *Vhist, *maskdbl, *u_out, *Vhist_out;
    int p, M, N, i, Tmax;
    int *mask;
    
    gamma = *mxGetPr(gamma_in);
    tol = *mxGetPr(tol_in);
    xtol = *mxGetPr(xtol_in);
    dt = *mxGetPr(dt_in);
    
    p = *mxGetPr(p_in);
    M = *mxGetPr(M_in);
    N = *mxGetPr(N_in);
    Tmax = *mxGetPr(Tmax_in);
    
    maskdbl = mxGetPr(mask_in);
    
    mask = (int*)mxCalloc(M*N, sizeof(int));
    for (i = 0; i < N*M; i++){
        mask[i] = (int)maskdbl[i];
    }
    
    g = mxGetPr(g_in);
    u = mxGetPr(u_in);
    dp = mxGetPr(dp_in);
    dqx = mxGetPr(dqx_in);
    dqy = mxGetPr(dqy_in);
    Suinvprod = mxGetPr(Suinvprod_in);
    Sginvprod = mxGetPr(Sginvprod_in); 
    Vhist = (double*)mxCalloc(Tmax+1,sizeof(double));
    
    residual = (double*)calloc(1,sizeof(double));
    
    // ------------------------------------------ evaluate function
    DGstep(g, u, mask, dp, dqx, dqy, Suinvprod, Sginvprod, p, gamma, tol, xtol, dt, M, N, Tmax, residual, Vhist);    
    
    // ------------------------------------------ deliver output [out, uT, dug, dxp, dxm, dyp, dym]
    const size_t dims[3] = {6, M, N};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    u_out = mxGetPr(plhs[0]);
    for (i = 0; i < 6*M*N; i++){
        u_out[i] = u[i];
    }
    plhs[1] = mxCreateDoubleScalar(*residual);
    plhs[2] = mxCreateDoubleMatrix(Tmax+1, 1, mxREAL);
    Vhist_out = mxGetPr(plhs[2]);
    for (i = 0; i < Tmax+1; i++){
        Vhist_out[i] = Vhist[i];
    }
    mxFree(Vhist);
}


void DGstep(double *g, double *u, int *mask, double *dp, double *dqx,
        double *dqy, double *Suinvprod, double *Sginvprod, int p, double gamma,
        double tol, double xtol, double dt, int M, int N, int Tmax, double *residual, double *Vhist){
    int i,j,k,l,m;
    double nrgtemp, nrgtemp2, Vo;
    double *out, *dug, *dxp, *dxm, *dyp, *dym;
    double gg[36];
    double uxm[36];
    double uym[36];
    double uxp[36];
    double uyp[36];
    int masks[4] = {0, 0, 0, 0};
    double x[2] = {-0.5, 0.5};
    double uT[6];
    double uTinv[6];
    double un[6];
    double V[9];
    double L[3];
    *residual = 1000;
    out = (double*)calloc(1,sizeof(double));
    dug = (double*)calloc(1,sizeof(double));
    dxp = (double*)calloc(1,sizeof(double));
    dxm = (double*)calloc(1,sizeof(double));
    dyp = (double*)calloc(1,sizeof(double));
    dym = (double*)calloc(1,sizeof(double));
    
    nrgtemp = 0;
    for (i = 0; i < M*N; i++){
        nrgtemp = nrgtemp + dp[i];
    }
    nrgtemp = nrgtemp/p;
    
    nrgtemp2 = 0;
    for (i = 0; i < M*(N+1); i++){
        nrgtemp2 = nrgtemp2 + dqx[i];
    }
    
    for (i = 0; i < (M+1)*N; i++){
        nrgtemp2 = nrgtemp2 + dqy[i];
    }
    nrgtemp2 *= gamma;
    Vhist[0] = nrgtemp + nrgtemp2;
    
    for (k = 0; k < Tmax; k++){
        if (*residual < tol){
            break;
        }
        for (i = 0; i < M; i++){
            for (j = 0; j < N; j++){
                if (mask[j + i*N]){
                    masks[0] = 0; masks[1] = 0; masks[2] = 0; masks[3] = 0;
                    if (i < M-1 && mask[j + (i+1)*N]){
                        masks[0] = 1;
                    }
                    if (i > 0 && mask[j + (i-1)*N]){
                        masks[1] = 1;
                    }
                    if (j < N-1 && mask[j+1 + N]){
                        masks[2] = 1;
                    }
                    if (j > 0 && mask[j-1 + i*N]){
                        masks[3] = 1;
                    }
                    
//                     printf("%i \n", j);
                    
                    for (l = 0; l < 6; l++){
                        uT[l] = u[6*j*M + 6*i + l];
                        un[l] = uT[l];
                    }
                    for (l = 0; l < 36; l++){
                        gg[l]  = Sginvprod[36*j*M + 36*i + l];
                        if (masks[0]){
                            uyp[l] = Suinvprod[36*j*M + 36*(i+1) + l];
                        }
                        if (masks[1]){
                            uym[l] = Suinvprod[36*j*M + 36*(i-1) + l];
                        }
                        if (masks[2]){
                            uxp[l] = Suinvprod[36*(j+1)*M + 36*i + l];
                        }
                        if (masks[3]){
                            uxm[l] = Suinvprod[36*(j-1)*M + 36*i + l];
                        }
                    }
                    
                    
                    for (l = 0; l < 6; l++){
                        inv3x3(uT,uTinv);
                        Vo = dp[i + j*M] + gamma*(dqx[i + j*M]+ dqx[i + (j+1)*M] + dqy[i + j*(M+1)] + dqy[i+1 + j*(M+1)]);
                        
                        *dug = 0; *dxp = 0; *dxm = 0; *dyp = 0; *dym = 0;
                        BrentDekker(x, uT, uTinv, uxp, uyp, uxm, uym, gg, p, gamma, l,
                                Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un, xtol, 1);
                        dp[i + j*M] = *dug;
                        
                        for (m = 0; m < 6; m++){
                            uT[m] = un[m];
                        }
//                         printf("%f \n", *dug);
                        if (masks[0]){
                            dqy[i+1 + j*(M+1)] = *dyp;
                        }
                        if (masks[1]){
                            dqy[i + j*(M+1)] = *dym;
                        }
                        if (masks[2]){
                            dqx[i + (j+1)*M] = *dxp;
                        }
                        if (masks[3]){
                            dqx[i + j*M] = *dxm;
                        }
                    }
                    
//                     for (m = 0; m < 6; m++){
//                         printf("%f \n", uT[m]);
//                     }
//                     printf("\n");

                    
                    eignvec3x3(uT,V,L);

                    for (l = 0; l < 6; l++){
                        u[6*j*M + 6*i + l] = uT[l];
                        
                    }
                    
                    prodform(V,L,gg);
                    
                    for (l = 0; l < 36; l++){
                        Suinvprod[36*j*M + 36*i + l] = gg[l];
                    }
                }
            }
        }
        
        
        nrgtemp = 0;
        for (i = 0; i < M*N; i++){
            nrgtemp = nrgtemp + dp[i];
        }
        nrgtemp = nrgtemp;
                
        nrgtemp2 = 0;
        for (i = 0; i < M*(N+1); i++){
            nrgtemp2 = nrgtemp2 + dqx[i];
        }
        
        for (i = 0; i < (M+1)*M; i++){
            nrgtemp2 = nrgtemp2 + dqy[i];
        }
        nrgtemp2 *= gamma;
        Vhist[k+1] = nrgtemp + nrgtemp2;
//         printf("%f \n", Vhist[k+1]);
        *residual = (Vhist[k] - Vhist[k+1])/Vhist[0];
    }

    
}

void eignvec3x3(double *A, double *V, double *L){
    double a1, a2, a3, a4, a5, a6, p1, p2, q, p, d, r, phi, eig1, eig2, eig3, norm1, norm2, norm3;
    int i,m;
    double A1[9];
    double A2[3];
    double A3[3];
    
    a1 = A[0];
    a2 = A[1];
    a3 = A[2];
    a4 = A[3];
    a5 = A[4];
    a6 = A[5];
    
    A1[0] = a1;
    A1[1] = a2;
    A1[2] = a3;
    A1[3] = a2;
    A1[4] = a4;
    A1[5] = a5;
    A1[6] = a3;
    A1[7] = a5;
    A1[8] = a6;
    
    p1 = a2*a2 + a3*a3 + a5*a5;
    q = (a1 + a4 + a6)/3;
    
    a1 = a1 - q;
    a4 = a4 - q;
    a6 = a6 - q;
    
    p2 = a1*a1 + a4*a4 + a6*a6 + 2*p1;
    p = sqrt(p2 / 6);
    
    d = a1*(a5*a5 - a4*a6) + a2*(a2*a6-a5*a3) + a3*(a4*a3 - a2*a5);
    
    r = -d/(2*p*p*p);
    
    if (r <= -1){
        phi = PI/3;
    }else if (r >= 1){
        phi = 0;
    }else{
        phi = acos(r)/3;
    }
            
    eig1 = q + 2*p*cos(phi);
    eig3 = q + 2*p*cos(phi + (2*PI/3));
    eig2 = 3*q - eig1 - eig3;
    
    for (i = 0; i < 3; i++){
        A2[i] = A1[i];
        A3[i] = A1[i];
    }
    
    A1[0] = A1[0] - eig2;
    A1[4] = A1[4] - eig2;
    A1[8] = A1[8] - eig2;
    
//     for (m = 0; m < 9; m++){
//         printf("%f \n", A1[m]);
//     }
    
    A2[0] = A2[0] - eig1;
    A3[0] = A3[0] - eig3;
    
    for (i = 0; i < 3; i++){
        V[i] = A1[i]*A3[0] + A1[i + 3]*A3[1] + A1[i + 6]*A3[2];
        V[i+6] = A1[i]*A2[0] + A1[i + 3]*A2[1] + A1[i + 6]*A2[2];
    }
    
    V[3] = V[1]*V[8] - V[2]*V[7];
    V[4] = V[2]*V[6] - V[0]*V[8];
    V[5] = V[0]*V[7] - V[1]*V[6];
    
    norm1 = 1.0/sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    norm2 = 1.0/sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
    norm3 = 1.0/sqrt(V[6]*V[6] + V[7]*V[7] + V[8]*V[8]);

    for (i = 0; i < 3; i++){
        V[i] = V[i]*norm1;
        V[i+3] = V[i+3]*norm2;
        V[i+6] = V[i+6]*norm3;
    }

    L[0] = 1/sqrt(eig1);
    L[1] = 1/sqrt(eig2);
    L[2] = 1/sqrt(eig3);
}

void prodform(double *V, double *L, double *S){
    double x1, x2, x3, x4, x5, x6, x1x2, x1x3, x2x2, x2x3, x2x4, x2x5, x3x3, x3x5, x3x6, x4x5, x5x5, x5x6, x1x5, x3x4, x2x6;
    double X[9];
    X[0] = L[0]*V[0]*V[0] + L[1]*V[3]*V[3] + L[2]*V[6]*V[6];
    X[1] = L[0]*V[0]*V[1] + L[1]*V[3]*V[4] + L[2]*V[6]*V[7];
    X[2] = L[0]*V[0]*V[2] + L[1]*V[3]*V[5] + L[2]*V[6]*V[8];
    X[3] = X[1];
    X[4] = L[0]*V[1]*V[1] + L[1]*V[4]*V[4] + L[2]*V[7]*V[7];
    X[5] = L[0]*V[1]*V[2] + L[1]*V[4]*V[5] + L[2]*V[7]*V[8];
    X[6] = X[2];
    X[7] = X[5];
    X[8] = L[0]*V[2]*V[2] + L[1]*V[5]*V[5] + L[2]*V[8]*V[8];
    
    x1 = X[0];
    x2 = X[1];
    x3 = X[2];
    x4 = X[4];
    x5 = X[5];
    x6 = X[8];
    
    x1x2 = x1*x2;
    x1x3 = x1*x3;
    
    x2x2 = x2*x2;
    x2x3 = x2*x3;
    x2x4 = x2*x4;
    x2x5 = x2*x5;
    
    x3x3 = x3*x3;
    x3x5 = x3*x5;
    x3x6 = x3*x6;
    
    x4x5 = x4*x5;
    
    x5x5 = x5*x5;
    x5x6 = x5*x6;
    
    x1x5 = x1*x5 + x2x3;
    x3x4 = x3*x4 + x2x5;
    x2x6 = x2*x6 + x3x5;
    
    S[0] = x1*x1; S[1] = 2*x1x2;       S[2] =  2*x1x3;       S[3] = x2x2;   S[4] = 2*x2x3;        S[5] = x3x3;
    S[6] = x1x2;  S[7] = x1*x4 + x2x2; S[8] =  x1x5;         S[9] = x2x4;   S[10] = x3x4;         S[11] = x3x5;
    S[12] = x1x3; S[13] = x1x5;        S[14] = x1*x6 + x3x3; S[15] = x2x5;  S[16] = x2x6;         S[17] = x3x6;
    S[18] = x2x2; S[19] = 2*x2x4;      S[20] = 2*x2x5;       S[21] = x4*x4; S[22] = 2*x4x5;       S[23] = x5x5;
    S[24] = x2x3; S[25] = x3x4;        S[26] = x2x6;         S[27] = x4x5;  S[28] = x4*x6 + x5x5; S[29] = x5x6;
    S[30] = x3x3; S[31] = 2*x3x5;      S[32] = 2*x3x6;       S[33] = x5x5;  S[34] = 2*x5x6;       S[35] = x6*x6;
    }

void inv3x3(double *uT, double* uTinv){
    double a1, a2, a3, a4, a5, a6, aa1, aa2, aa3, d;
    a1 = uT[0];
    a2 = uT[1];
    a3 = uT[2];
    a4 = uT[3];
    a5 = uT[4];
    a6 = uT[5];
    aa1 = a6*a4-a5*a5;
    aa2 = a5*a3-a2*a6;
    aa3 = a5*a2-a4*a3;
    
    d = 1/(a1*aa1 + a2*aa2 + a3*aa3);
    
    uTinv[0] = d*aa1;
    uTinv[1] = d*aa2;
    uTinv[2] = d*aa3;
    uTinv[3] = d*(a6*a1 - a3*a3);
    uTinv[4] = d*(a3*a2 - a5*a1);
    uTinv[5] = d*(a4*a1 - a2*a2);
    }

void BrentDekker(double *x, double *uT, double *uTinv, double *uxp, double *uyp,
        double *uxm, double *uym, double *g, int pp, double gamma, int l, double Vo, double dt, int *masks,
        double *out, double* dug, double *dxp, double *dxm, double* dyp, double* dym, double* un, double tol, int interval){
    bool badint;
    double a, b, fa, fb, fx, dx, twosqrt,fc, c, d, e, m, p, q, r, s, toler, xx;
    int i;
    badint = 0;
            xx = 0.01;
    if (interval) {
//         for (i = 0; i < 36; i++){
//             printf("%f \n",g[i]);
//        }
        a = x[0];
        b = x[1];
        Vdiff_solv(a, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
//         printf("%f \n", *dug);
//         printf("%f \n", *dxp);
//         printf("%f \n", *dxm);
//         printf("%f \n", *dyp);
//         printf("%f \n", *dym);
//         printf("\n");
        fa = *out;
        Vdiff_solv(b, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
        fb = *out;

        
        if (fa*fb > 0){
            badint = 1;
        }
        
        if (!fa){
            b = a;
            return;
        }else if (!fb){
            return;
        }
    }
    if (!interval || badint){
        Vdiff_solv(xx, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
        fx = *out;
        
        if (fx == 0){
            b = xx;
            return;
        }
        if (x != 0){
            dx = xx/50;
        } else {
            dx = 1/50;
        }
        
        twosqrt = sqrt(2);
        a = xx; fa = fx; b = xx; fb = fx;
        
        while ((fa > 0) == (fb > 0)){
            dx = twosqrt*dx;
            a = xx - dx;  
            Vdiff_solv(a, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
            fa = *out;
            if ((fa > 0) != (fb > 0)){
                break;
            }
            
            b = xx + dx;
            Vdiff_solv(b, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
            fb = *out;
        }
    }
    fc = fb;

    while (fb != 0 && a != b){
        if ((fb > 0) == (fc > 0)){
            c = a;  fc = fa;
            d = b - a;  e = d;
        }
        if (fabs(fc) < fabs(fb)){
            a = b;    b = c;    c = a;
            fa = fb;  fb = fc;  fc = fa;
            *out = fb;
        }
        
        m = 0.5*(c - b);
        if (fabs(b) > 1.0){
            toler = 2.0*tol*fabs(b);
        }else{
            toler = 2.0*tol;
        }
        if ((fabs(m) <= toler) || (fb == 0.0)){
            break;
        }
        
        
        if ((fabs(e) < toler) || (fabs(fa) <= fabs(fb))){
            d = m;  e = m;
        } else{
            s = fb/fa;
            if (a == c){
                p = 2.0*m*s;
                q = 1.0 - s;
            } else{
                q = fa/fc;
                r = fb/fc;
                p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
                q = (q - 1.0)*(r - 1.0)*(s - 1.0);
            }
            if (p > 0){
                q = -q;
            } else{
                p = -p;
            }
            if ((2.0*p < 3.0*m*q - fabs(toler*q)) && (p < fabs(0.5*e*q))){
                e = d;  d = p/q;
            }else{
                d = m;  e = m;
            }
        }
        a = b;
        fa = fb;
        if (fabs(d) > toler){
            b = b + d;
        }else if (b > c){
            b = b - toler;
        }else {
            b = b + toler;
        }
        
        Vdiff_solv(b, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
        fb = *out;
//         printf("%e \n",fb);
    }
}



void Vdiff_solv(double alpha, double *uT, double *uTinv, double *uxp, double *uyp,
        double *uxm, double *uym, double *g, int p, double gamma, int l, double Vo, double dt, int *masks,
        double *out, double* dug, double *dxp, double *dxm, double* dyp, double* dym, double* un){
    
    // ------------------------------------------ initialization
    double C;
    int i;
    
    C = 0.5*alpha*alpha;
    if (l == 0){
        un[0] = uT[0] + alpha + C*uTinv[0];
    }else if (l == 1){
        un[3] = uT[3] + alpha + C*uTinv[3];
    }else if (l == 2){
        un[5] = uT[5] + alpha + C*uTinv[5];
    }else if (l == 3){
        un[0] = uT[0] + C*uTinv[3];
        un[1] = uT[1] + alpha + C*uTinv[1];
        un[3] = uT[3] + C*uTinv[0];
    }else if (l == 4){
        un[3] = uT[3] + C*uTinv[5];
        un[4] = uT[4] + alpha + C*uTinv[4];
        un[5] = uT[5] + C*uTinv[3];
    }else{
        un[0] = uT[0] + C*uTinv[5];
        un[2] = uT[2] + alpha + C*uTinv[2];
        un[5] = uT[5] + C*uTinv[0];
    }
    
//     for (i = 0; i < 6; i++){
//         printf("%f \n",un[i]);
//     }
//     printf("\n");
            
    if (p == 1){
        *dug = d_M(g,un);
        
    }else{
        *dug = d_M(g,un);
        *dug = (*dug)*(*dug)/p;
    }
        
    if (masks[0]){
        *dyp = d_M(uyp,un);
    }
    if (masks[1]){
        *dym = d_M(uym,un);
    }
    if (masks[2]){
        *dxp = d_M(uxp,un);
    }
    if (masks[3]){
        *dxm = d_M(uxm,un);
    }
            
    *out = *dug - Vo + gamma*(*dxp + *dxm + *dyp + *dym);
    *out = alpha + dt*(*out)/alpha;
    
}

double d_M(double *XX, double *A){
    double A1, A2, A3, A4, A5, A6, p1, q, p2, p, d, r, eig1, eig2, eig3, phi;

    A1 = XX[0]*A[0]  + XX[1]*A[1]  + XX[2]*A[2]  + XX[3]*A[3]  + XX[4]*A[4] + XX[5]*A[5];
    A2 = XX[6]*A[0]  + XX[7]*A[1]  + XX[8]*A[2]  + XX[9]*A[3]  + XX[10]*A[4] + XX[11]*A[5];
    A3 = XX[12]*A[0] + XX[13]*A[1] + XX[14]*A[2] + XX[15]*A[3] + XX[16]*A[4] + XX[17]*A[5];
    A4 = XX[18]*A[0] + XX[19]*A[1] + XX[20]*A[2] + XX[21]*A[3] + XX[22]*A[4] + XX[23]*A[5];
    A5 = XX[24]*A[0] + XX[25]*A[1] + XX[26]*A[2] + XX[27]*A[3] + XX[28]*A[4] + XX[29]*A[5];
    A6 = XX[30]*A[0] + XX[31]*A[1] + XX[32]*A[2] + XX[33]*A[3] + XX[34]*A[4] + XX[35]*A[5];
    
    p1 = A2*A2 + A3*A3 + A5*A5;
    q = (A1 + A4 + A6)/3;
    
    A1 = A1 - q;
    A4 = A4 - q;
    A6 = A6 - q;
    
    p2 = A1*A1 + A4*A4 + A6*A6 + 2*p1;
    p = 2*sqrt(p2 / 6);
    
    d = A5*(A1*A5 - 2*A2*A3) + A6*(A2*A2 - A1*A4) + A4*A3*A3;
    r = -4*d/(p*p*p);
    
    if (r < -1){
        eig1 = q + 0.5*p;
        eig3 = q - p;
    }else if (r > 1){
        eig1 = q + p;
        eig3 = q - 0.5*p;
    }else{
        phi = acos(r)/3;
        eig1 = q + p*cos(phi);
        eig3 = q + p*cos(phi + 2*PI/3);
    }
    
    eig2 = 3*q - eig1 - eig3;
    
    eig1 = log(eig1);
    eig2 = log(eig2);
    eig3 = log(eig3);
    
    return eig1*eig1 + eig2*eig2 + eig3*eig3;
    
}
