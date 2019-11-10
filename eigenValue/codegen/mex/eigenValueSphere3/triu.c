/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * triu.c
 *
 * Code generation for function 'triu'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "eigenValueSphere3.h"
#include "triu.h"

/* Function Definitions */
void triu(real_T x[160000])
{
  int32_T istart;
  int32_T j;
  int32_T i;
  istart = 1;
  for (j = 0; j < 400; j++) {
    for (i = istart; i < 401; i++) {
      x[(i + 400 * j) - 1] = 0.0;
    }

    istart++;
  }
}

/* End of code generation (triu.c) */
