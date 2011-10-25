#include "getEs.h" 
#include <R.h>
#include <Rinternals.h>

double absolute(double x)
{
  if (x<0)
  return (-x);
  else
  return (x);
}

void cumsum(double *x, int len)
{
  int i;
  for (i = 1; i < len; ++i) {
    *(x + i) = *(x + i) + *(x + i -1);
  }
}

double getNr(double *fchr, int *sign, int signLen)
{
  int i;
  double nr;
  nr = 0.0;
  for (i = 0; i < signLen; ++i) {
    nr = absolute(fchr[sign[i] -1]) + nr;
    }
  return nr;
}

void getPhit(double *fchr, int *sign, int signLen, double nr, double *phit)
{
  int i;
  for (i = 0; i < signLen; ++i) {
    *(phit + sign[i]-1) = absolute(*(fchr + sign[i]-1)) / nr;
    }
}

void getPmiss(int *sign, int fchrLen, int signLen, double *pmiss)
{
  int i;
  double tmp = 1.0 / (fchrLen-signLen);
  for (i = 0; i < fchrLen; ++i) {
    *(pmiss + i) = tmp;
    }
  for (i = 0; i < signLen; ++i) {
    *(pmiss + sign[i]-1) = 0;
    }
}

SEXP getEs(SEXP fchr, SEXP sign)
{ 
  PROTECT(fchr = coerceVector(fchr, REALSXP)); 
  PROTECT(sign = coerceVector(sign, INTSXP));

  double *rfchr = REAL(fchr);
  int *rsign = INTEGER(sign);

  int i, nfchr, nsign;

  nfchr = length(fchr);
  nsign = length(sign);

  SEXP es;
  PROTECT(es = allocVector(REALSXP, nfchr));
  double *res;
  res = REAL(es);
  for(i = 0; i < nfchr; i++) res[i] = 0.0;

  double nr = getNr(rfchr, rsign, nsign);

  SEXP phit;
  PROTECT(phit = allocVector(REALSXP, nfchr));
  double *rphit;
  rphit = REAL(phit);
  for(i = 0; i < nfchr; i++) rphit[i] = 0.0;
  getPhit(rfchr, rsign, nsign, nr, rphit);
  cumsum(rphit, nfchr);

  SEXP pmiss;
  PROTECT(pmiss = allocVector(REALSXP, nfchr));
  double *rpmiss;
  rpmiss = REAL(pmiss);
  for(i = 0; i < nfchr; i++) rpmiss[i] = 0.0;
  getPmiss(rsign, nfchr, nsign, rpmiss);
  cumsum(rpmiss, nfchr);

  for (i = 0; i < nfchr; ++i) {
    res[i] = rphit[i] - rpmiss[i];
    }

  UNPROTECT(5);
  return es;
}
