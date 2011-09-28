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
  for (i = 0; i < len; ++i) {
    *(x + i) = *(x + i) + *(x + i -1);
  }
}

double getNr(double *fchr, int *sign, int signLen)
{
  int i;
  double nr;
  nr = 0;
  for (i = 0; i < signLen; ++i) {
    nr = absolute(fchr[sign[i]-1]) + nr;
    }
  return nr;
}

void populateArray(double *xin, int len)
{
  int i;
  for (i = 0; i < len; ++i) {
    *(xin + i) = 0.0;
    }
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
  int i, nfchr, nsign;
  double *rfchr = REAL(fchr), *res;
  int *rsign = INTEGER(sign);

  PROTECT(fchr = coerceVector(fchr, REALSXP)); 
  PROTECT(sign = coerceVector(sign, REALSXP));

  nfchr = length(fchr);
  nsign = length(sign);

  SEXP es;
  PROTECT(es = allocVector(REALSXP, nfchr));
  res = REAL(es);

  double nr = getNr(rfchr, rsign, nsign);

  double phit[nfchr];
  double *rphit = calloc(nfchr, sizeof(double));
  populateArray(phit,nfchr);
  getPhit(rfchr, rsign, nsign, nr, phit);
  cumsum(phit, nfchr);

  double pmiss[nfchr];
  double *rpmiss = calloc(nfchr, sizeof(double));
  populateArray(pmiss,nfchr);
  getPmiss(rsign, nfchr, nsign, pmiss);
  cumsum(pmiss, nfchr);

  for (i = 0; i < nfchr; ++i) {
    res[i] = phit[i] - pmiss[i];
    }

  UNPROTECT(3);
  return es;
}
