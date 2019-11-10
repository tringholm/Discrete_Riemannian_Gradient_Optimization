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

#define x_in             prhs[0]
#define uT_in            prhs[1]
#define uTinv_in         prhs[2]
#define uxp_in           prhs[3]
#define uyp_in           prhs[4]
#define uxm_in           prhs[5]
#define uym_in           prhs[6]
#define g_in             prhs[7]
#define p_in             prhs[8]
#define gamma_in         prhs[9]
#define l_in             prhs[10]
#define Vo_in            prhs[11]
#define dt_in            prhs[12]
#define masks_in         prhs[13]
#define tol_in           prhs[14]
#define interval_in      prhs[15]

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
    double gamma, Vo, dt, tol;
    double *uT, *uTinv, *uxp, *uyp, *uxm, *uym, *g, *out, *dug, *dxp, *dxm, *dyp, *dym, *u_out, *masksdbl, *x;
    int p, l, i, interval;
    int masks[4];
    double un[6];
    
    x = mxGetPr(x_in);
    gamma = *mxGetPr(gamma_in);
    Vo = *mxGetPr(Vo_in);
    dt = *mxGetPr(dt_in);
    tol = *mxGetPr(tol_in);
    
    p = *mxGetPr(p_in);
    l = *mxGetPr(l_in);
    interval = *mxGetPr(interval_in);
    
    masksdbl = mxGetPr(masks_in);
    masks[0] = (int)masksdbl[0];
    masks[1] = (int)masksdbl[1];
    masks[2] = (int)masksdbl[2];
    masks[3] = (int)masksdbl[3];
    
    uT = mxGetPr(uT_in);
    uTinv = mxGetPr(uTinv_in);
    uxp = mxGetPr(uxp_in);
    uyp = mxGetPr(uyp_in);
    uxm = mxGetPr(uxm_in);
    uym = mxGetPr(uym_in);
    g = mxGetPr(g_in); 
    
    for (i = 0; i < 6; i++){
        un[i] = uT[i];
    }
    
    out = (double*)calloc(1,sizeof(double));
    dug = (double*)calloc(1,sizeof(double));
    dxp = (double*)calloc(1,sizeof(double));
    dxm = (double*)calloc(1,sizeof(double));
    dyp = (double*)calloc(1,sizeof(double));
    dym = (double*)calloc(1,sizeof(double));
    
    // ------------------------------------------ evaluate function
    BrentDekker(x, uT, uTinv, uxp, uyp, uxm, uym, g, p, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un,tol,interval);
    
    
    // ------------------------------------------ deliver output [out, uT, dug, dxp, dxm, dyp, dym]
    plhs[0] = mxCreateDoubleScalar(*out);
    plhs[1] = mxCreateDoubleMatrix(6, 1,mxREAL);
    u_out = mxGetPr(plhs[1]);
    for (i = 0; i < 6; i++){
        u_out[i] = un[i];
    }
    plhs[2] = mxCreateDoubleScalar(*dug);
    plhs[3] = mxCreateDoubleScalar(*dxp);
    plhs[4] = mxCreateDoubleScalar(*dxm);
    plhs[5] = mxCreateDoubleScalar(*dyp);
    plhs[6] = mxCreateDoubleScalar(*dym);
}


void BrentDekker(double *x, double *uT, double *uTinv, double *uxp, double *uyp,
        double *uxm, double *uym, double *g, int pp, double gamma, int l, double Vo, double dt, int *masks,
        double *out, double* dug, double *dxp, double *dxm, double* dyp, double* dym, double* un, double tol, int interval){
    bool badint;
    double a, b, fa, fb, fx, dx, twosqrt,fc, c, d, e, m, p, q, r, s, toler, xx;
    badint = 0;
            xx = 0.01;
    if (interval) {
        
        a = x[0];
        b = x[1];
        Vdiff_solv(a, uT, uTinv, uxp, uyp, uxm, uym, g, pp, gamma, l, Vo, dt, masks, out, dug, dxp, dxm, dyp, dym, un);
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
    
    C = 0.5*alpha*alpha;
    if (l == 1){
        un[0] = uT[0] + alpha + C*uTinv[0];
    }else if (l == 2){
        un[3] = uT[3] + alpha + C*uTinv[3];
    }else if (l == 3){
        un[5] = uT[5] + alpha + C*uTinv[5];
    }else if (l == 4){
        un[0] = uT[0] + C*uTinv[3];
        un[1] = uT[1] + alpha + C*uTinv[1];
        un[3] = uT[3] + C*uTinv[0];
    }else if (l == 5){
        un[3] = uT[3] + C*uTinv[5];
        un[4] = uT[4] + alpha + C*uTinv[4];
        un[5] = uT[5] + C*uTinv[3];
    }else{
        un[0] = uT[0] + C*uTinv[5];
        un[2] = uT[2] + alpha + C*uTinv[2];
        un[5] = uT[5] + C*uTinv[0];
    }
            
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
