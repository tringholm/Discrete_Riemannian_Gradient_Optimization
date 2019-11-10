/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eigenValueSphere3_terminate.c
 *
 * Code generation for function 'eigenValueSphere3_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "eigenValueSphere3.h"
#include "eigenValueSphere3_terminate.h"
#include "_coder_eigenValueSphere3_mex.h"
#include "eigenValueSphere3_data.h"

/* Function Definitions */
void eigenValueSphere3_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void eigenValueSphere3_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (eigenValueSphere3_terminate.c) */
