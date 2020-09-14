/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: succLam.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 28-Feb-2020 15:32:32
 */

/* Include Files */
#include "succLam.h"
#include "succLam_emxutil.h"

/* Function Declarations */
static int div_s32_floor(int numerator, int denominator);

/* Function Definitions */

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_s32_floor(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  boolean_T quotientNeedsNegation;
  unsigned int tempAbsQuotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~(unsigned int)numerator + 1U;
    } else {
      absNumerator = (unsigned int)numerator;
    }

    if (denominator < 0) {
      absDenominator = ~(unsigned int)denominator + 1U;
    } else {
      absDenominator = (unsigned int)denominator;
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }

      quotient = -(int)tempAbsQuotient;
    } else {
      quotient = (int)tempAbsQuotient;
    }
  }

  return quotient;
}

/*
 * Compute successive total infectiousness for this I
 * Arguments    : const emxArray_real_T *Icurr
 *                const emxArray_real_T *Pomega
 *                double ncurr
 *                emxArray_real_T *Lsucc
 * Return Type  : void
 */
void succLam(const emxArray_real_T *Icurr, const emxArray_real_T *Pomega, double
             ncurr, emxArray_real_T *Lsucc)
{
  int i;
  int loop_ub;
  int b_i;
  double b_Icurr;
  int i1;

  /*  Simple function to compute multiplication in turn */
  i = Lsucc->size[0] * Lsucc->size[1];
  Lsucc->size[0] = 1;
  loop_ub = (int)ncurr;
  Lsucc->size[1] = loop_ub;
  emxEnsureCapacity_real_T(Lsucc, i);
  for (i = 0; i < loop_ub; i++) {
    Lsucc->data[i] = 0.0;
  }

  i = (int)(ncurr + -1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /*  Relevant part of SI: Pomega(1:i-1)) */
    loop_ub = div_s32_floor(-b_i, -1);
    if ((loop_ub + 1 == 1) || (b_i + 1 == 1)) {
      b_Icurr = 0.0;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_Icurr += Icurr->data[b_i - i1] * Pomega->data[i1];
      }

      Lsucc->data[b_i + 1] = b_Icurr;
    } else {
      b_Icurr = 0.0;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_Icurr += Icurr->data[b_i - i1] * Pomega->data[i1];
      }

      Lsucc->data[b_i + 1] = b_Icurr;
    }
  }
}

/*
 * File trailer for succLam.c
 *
 * [EOF]
 */
