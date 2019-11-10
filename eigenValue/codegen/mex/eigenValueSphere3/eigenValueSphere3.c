/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eigenValueSphere3.c
 *
 * Code generation for function 'eigenValueSphere3'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "eigenValueSphere3.h"
#include "eigenValueSphere3_emxutil.h"
#include "mpower.h"
#include "sum.h"
#include "triu.h"
#include "diag.h"
#include "blas.h"

/* Function Declarations */
static void euclidCoords(const real_T phi[399], real_T x[400]);

/* Function Definitions */
static void euclidCoords(const real_T phi[399], real_T x[400])
{
  real_T sinefact;
  int32_T i;
  x[0] = muDoubleScalarCos(phi[0]);
  sinefact = muDoubleScalarSin(phi[0]);
  for (i = 0; i < 398; i++) {
    x[i + 1] = muDoubleScalarCos(phi[i + 1]) * sinefact;
    sinefact *= muDoubleScalarSin(phi[i + 1]);
  }

  x[399] = sinefact;
}

void eigenValueSphere3(eigenValueSphere3StackData *SD, const real_T A[160000],
  real_T phi[399], real_T tol, real_T dt, emxArray_real_T *Vhist)
{
  real_T residual;
  real_T dA[400];
  real_T b_Vhist[2000];
  real_T u[400];
  real_T cpa;
  real_T Vd;
  char_T TRANSB;
  char_T TRANSA;
  real_T rowsums[400];
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  int32_T cnt;
  real_T sums_idx_0;
  real_T sums_idx_1;
  real_T sums_idx_2;
  real_T sums_idx_3;
  real_T sums_idx_4;
  int32_T i;
  real_T sum3fact;
  int32_T j;
  real_T ufact;
  real_T alpha;
  real_T dta;
  real_T cp;
  real_T sp;
  real_T cc;
  real_T ss;
  real_T sc;
  real_T V;
  real_T b_u[400];

  /*  Compute the lowest eigenvalue/vector of a symmetric matrix by minimizing the */
  /*  Rayleigh quotient over S^{n-1}. */
  residual = rtInf;
  diag(A, dA);
  memcpy(&SD->f0.AU[0], &A[0], 160000U * sizeof(real_T));
  triu(SD->f0.AU);
  memset(&b_Vhist[0], 0, 2000U * sizeof(real_T));
  euclidCoords(phi, u);
  cpa = 1.0;
  Vd = 0.0;
  TRANSB = 'N';
  TRANSA = 'N';
  memset(&rowsums[0], 0, 400U * sizeof(real_T));
  m_t = (ptrdiff_t)400;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)400;
  lda_t = (ptrdiff_t)400;
  ldb_t = (ptrdiff_t)400;
  ldc_t = (ptrdiff_t)400;
  dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &cpa, &A[0], &lda_t, &u[0], &ldb_t,
        &Vd, &rowsums[0], &ldc_t);
  for (cnt = 0; cnt < 400; cnt++) {
    rowsums[cnt] *= u[cnt];
  }

  b_Vhist[0] = sum(rowsums);

  /*  tic */
  sums_idx_0 = A[0] * mpower(u[0]);

  /*  i = k, j = k */
  sums_idx_1 = 0.0;

  /*  i < k, j = k */
  sums_idx_2 = 2.0 * (rowsums[0] - sums_idx_0);

  /*  i = k, j > k */
  sums_idx_3 = 0.0;

  /*  i < k, j > k */
  sums_idx_4 = (b_Vhist[0] - sums_idx_2) - sums_idx_0;

  /*  i > k, j > k */
  cpa = 1.0;
  Vd = 0.0;
  TRANSB = 'N';
  TRANSA = 'N';
  memset(&rowsums[0], 0, 400U * sizeof(real_T));
  m_t = (ptrdiff_t)400;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)400;
  lda_t = (ptrdiff_t)400;
  ldb_t = (ptrdiff_t)400;
  ldc_t = (ptrdiff_t)400;
  dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &cpa, &SD->f0.AU[0], &lda_t, &u[0],
        &ldb_t, &Vd, &rowsums[0], &ldc_t);
  for (cnt = 0; cnt < 400; cnt++) {
    rowsums[cnt] *= 2.0 * u[cnt];
  }

  i = 1;
  while ((i - 1 < 2000) && (!(residual < tol))) {
    sum3fact = 1.0;
    ufact = 1.0;
    for (j = 0; j < 399; j++) {
      alpha = 0.1;
      dta = dt / 0.1;
      cp = 1.0 / muDoubleScalarCos(phi[j]);
      sp = 1.0 / muDoubleScalarSin(phi[j]);
      cpa = muDoubleScalarCos(phi[j] + 0.1);
      residual = muDoubleScalarSin(phi[j] + 0.1);
      cc = cpa * cp;
      ss = residual * sp;
      sc = residual * cp;
      residual = cpa * sp;
      Vd = (((sums_idx_0 * (mpower(cc) - 1.0) + sums_idx_1 * (cc - 1.0)) +
             sums_idx_2 * (cc * ss - 1.0)) + sums_idx_3 * (ss - 1.0)) +
        sums_idx_4 * (mpower(ss) - 1.0);
      V = 0.1 + dta * Vd;
      cpa = rtInf;
      cnt = 0;
      while ((muDoubleScalarAbs(V) > 1.0E-13) && (cnt < 10) && (cpa - V > 0.0))
      {
        alpha -= V / ((1.0 + dta * ((((sums_idx_0 * (-2.0 * sc * cc) +
          sums_idx_1 * -sc) + sums_idx_2 * (cc * residual - sc * ss)) +
          sums_idx_3 * residual) + sums_idx_4 * 2.0 * residual * ss)) - dta * Vd
                      / alpha);
        dta = dt / alpha;
        cpa = muDoubleScalarCos(phi[j] + alpha);
        residual = muDoubleScalarSin(phi[j] + alpha);
        cc = cpa * cp;
        ss = residual * sp;
        sc = residual * cp;
        residual = cpa * sp;
        Vd = (((sums_idx_0 * (mpower(cc) - 1.0) + sums_idx_1 * (cc - 1.0)) +
               sums_idx_2 * (cc * ss - 1.0)) + sums_idx_3 * (ss - 1.0)) +
          sums_idx_4 * (mpower(ss) - 1.0);
        cpa = V;
        V = alpha + dta * Vd;
        cnt++;
      }

      u[j] = u[j] * cc * ufact;
      ufact *= ss;
      residual = u[j + 1] * ufact;
      sum3fact *= mpower(ss);
      sums_idx_0 = residual * A[(j + 400 * (j + 1)) + 1] * residual;

      /*  i = k, j = k */
      sums_idx_1 = 0.0;
      for (cnt = 0; cnt <= j; cnt++) {
        sums_idx_1 += A[cnt + 400 * (j + 1)] * u[cnt];
      }

      sums_idx_1 *= 2.0 * residual;
      sums_idx_3 = (cc * ss * sums_idx_2 + ss * sums_idx_3) - sums_idx_1;

      /*  i < k, j > k */
      sums_idx_2 = rowsums[j + 1] * sum3fact;
      sums_idx_4 = (ss * ss * sums_idx_4 - sums_idx_2) - sums_idx_0;

      /*  i > k, j > k */
      phi[j] += alpha;
    }

    u[399] *= ufact;
    cpa = 1.0;
    Vd = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    memset(&rowsums[0], 0, 400U * sizeof(real_T));
    m_t = (ptrdiff_t)400;
    n_t = (ptrdiff_t)1;
    k_t = (ptrdiff_t)400;
    lda_t = (ptrdiff_t)400;
    ldb_t = (ptrdiff_t)400;
    ldc_t = (ptrdiff_t)400;
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &cpa, &SD->f0.AU[0], &lda_t, &u[0],
          &ldb_t, &Vd, &rowsums[0], &ldc_t);
    for (cnt = 0; cnt < 400; cnt++) {
      b_u[cnt] = u[cnt] * dA[cnt] * u[cnt];
      rowsums[cnt] *= 2.0 * u[cnt];
    }

    b_Vhist[i] = sum(rowsums) + sum(b_u);
    sums_idx_0 = A[0] * mpower(u[0]);

    /*  i = k, j = k */
    sums_idx_1 = 0.0;

    /*  i < k, j = k */
    sums_idx_2 = rowsums[0];

    /*  i = k, j > k */
    sums_idx_3 = 0.0;

    /*  i < k, j > k */
    sums_idx_4 = (b_Vhist[i] - rowsums[0]) - sums_idx_0;

    /*  i > k, j > k */
    residual = (b_Vhist[i - 1] - b_Vhist[i]) / muDoubleScalarAbs(b_Vhist[0]);
    i++;
  }

  /*  toc */
  j = 0;
  for (i = 0; i < 2000; i++) {
    if (b_Vhist[i] != 0.0) {
      j++;
    }
  }

  cnt = Vhist->size[0];
  Vhist->size[0] = j;
  emxEnsureCapacity((emxArray__common *)Vhist, cnt, sizeof(real_T));
  cnt = 0;
  for (i = 0; i < 2000; i++) {
    if (b_Vhist[i] != 0.0) {
      Vhist->data[cnt] = b_Vhist[i];
      cnt++;
    }
  }
}

/* End of code generation (eigenValueSphere3.c) */
