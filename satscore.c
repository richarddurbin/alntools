/*  File: satscore.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: scoring shared by satscan and satmatch, see satscore.h
 * Exported functions: see satscore.h
 * HISTORY:
 * Last edited: Oct 26 17:43 2026 (claude)
 * Created: Mon Oct 26 12:00:00 2026 (claude)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "satscore.h"

const SatScheme satScheme[MAX_SCHEME] = {
  { 2, -1, -3, -15, -1 },	// weak
  { 2, -2, -4, -20, -2 },	// medium
  { 2, -3, -6, -25, -2 } } ;	// strong
const char satSchemeChar[MAX_SCHEME] = { 'w', 'm', 's' } ;

double satExpectedScore (const SatScheme *sp, const double *p, const double *q)
{
  double e = 0 ;
  for (int a = 0 ; a < 4 ; ++a)
    for (int b = 0 ; b < 4 ; ++b) e += p[a] * q[b] * satBaseScore (sp, a, b) ;
  return e ;
}

int satSchemeLevel (const char *seq, int len)
{
  double p[4] = { 0, 0, 0, 0 } ;
  for (int i = 0 ; i < len ; ++i) p[(int)seq[i]] += 1.0 / len ;
  int s ;
  for (s = 0 ; s < MAX_SCHEME ; ++s)
    if (satExpectedScore (&satScheme[s], p, p) <= DANGER_T) break ;
  return s ;
}

// Rows are the bases of t, columns the bases of q. h[j] is the best score of an alignment of t so far
// whose last aligned base of q is q[j]; f[j] the best ending in a gap in q at j (bases of t inserted);
// within a row a gap in t (bases of q deleted) runs along the row and can wrap from q[qLen-1] to q[0].
// Before the first row h is 0 everywhere, so the alignment can start at any rotation of q.
int satGlobalScore (const char *t, int tLen, const char *q, int qLen, const SatScheme *sp, int threshold)
{
  static const int MINUS_INF = -(1 << 28) ;
  int *prev = new (qLen, int), *h = new (qLen, int), *f = new (qLen, int), *gIn = new (qLen, int) ;
  for (int j = 0 ; j < qLen ; ++j) { prev[j] = 0 ; f[j] = MINUS_INF ; }
  int best = MINUS_INF ;
  for (int i = 0 ; i < tLen ; ++i)
    { int c = t[i] ;
      for (int j = 0 ; j < qLen ; ++j) // diagonal from the previous row, and gaps in q from above
	{ int d = prev[j ? j-1 : qLen-1] + satBaseScore (sp, c, q[j]) ;
	  int o = prev[j] + sp->gapOpen, e = f[j] + sp->gapExtend ;
	  f[j] = o > e ? o : e ;
	  h[j] = d > f[j] ? d : f[j] ;
	}
      // gaps in t along the row: a first pass from q[0], entering with a gap opened from q[qLen-1], then
      // continuing round the wrap while the gap beats the one that arrived at each cell on the first pass
      int g = h[qLen-1] + sp->gapOpen ;
      for (int j = 0 ; j < qLen ; ++j)
	{ gIn[j] = g ;
	  if (g > h[j]) { h[j] = g ; g += sp->gapExtend ; }
	  else { int o = h[j] + sp->gapOpen, e = g + sp->gapExtend ; g = o > e ? o : e ; }
	}
      for (int j = 0 ; j < qLen && g > gIn[j] ; ++j, g += sp->gapExtend)
	if (g > h[j]) h[j] = g ;
      int rowMax = MINUS_INF ;
      for (int j = 0 ; j < qLen ; ++j) if (h[j] > rowMax) rowMax = h[j] ;
      int *tmp = prev ; prev = h ; h = tmp ;
      best = rowMax ;
      if (rowMax + sp->match * (tLen - 1 - i) <= threshold) // can no longer end above threshold
	{ best = rowMax + sp->match * (tLen - 1 - i) ; break ; }
    }
  newFree (prev, qLen, int) ; newFree (h, qLen, int) ; newFree (f, qLen, int) ; newFree (gIn, qLen, int) ;
  return best ;
}

/******************* end of file **************/
