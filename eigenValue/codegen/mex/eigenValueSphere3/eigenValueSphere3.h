/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eigenValueSphere3.h
 *
 * Code generation for function 'eigenValueSphere3'
 *
 */

#ifndef EIGENVALUESPHERE3_H
#define EIGENVALUESPHERE3_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "eigenValueSphere3_types.h"

/* Function Declarations */
extern void eigenValueSphere3(eigenValueSphere3StackData *SD, const real_T A
  [160000], real_T phi[399], real_T tol, real_T dt, emxArray_real_T *Vhist);

#endif

/* End of code generation (eigenValueSphere3.h) */
