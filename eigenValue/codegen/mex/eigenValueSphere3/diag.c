/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diag.c
 *
 * Code generation for function 'diag'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "eigenValueSphere3.h"
#include "diag.h"

/* Function Definitions */
void diag(const real_T v[160000], real_T d[400])
{
  int32_T j;
  for (j = 0; j < 400; j++) {
    d[j] = v[j * 401];
  }
}

/* End of code generation (diag.c) */
