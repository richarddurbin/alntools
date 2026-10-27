/*  File: satscore.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: scoring shared by satscan and satmatch, see satscore.h
 * Exported functions: see satscore.h
 * HISTORY:
 * Last edited: Oct 27 01:37 2026 (claude)
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

// an ordinary global affine alignment of all of s against all of r, starting at column start of r
static int rotatedScore (const char *s, int n, const char *r, int L, int start, const SatScheme *sp, int threshold)
{
  static const int MINUS_INF = -(1 << 28) ;
  int *h = new (L+1, int), *hNew = new (L+1, int), *f = new (L+1, int) ;
  h[0] = 0 ; f[0] = MINUS_INF ;
  for (int j = 1 ; j <= L ; ++j) { h[j] = sp->gapOpen + (j-1) * sp->gapExtend ; f[j] = MINUS_INF ; }
  int result = h[L] ;
  for (int i = 0 ; i < n ; ++i)
    { hNew[0] = sp->gapOpen + i * sp->gapExtend ;
      int e = MINUS_INF, rowMax = hNew[0] ;
      for (int j = 1, k = start ; j <= L ; ++j, k = (k+1 == L) ? 0 : k+1)
	{ int o = h[j] + sp->gapOpen, x = f[j] + sp->gapExtend ; f[j] = o > x ? o : x ;	// s[i] alone
	  o = hNew[j-1] + sp->gapOpen ; x = e + sp->gapExtend ; e = o > x ? o : x ;	// r[k] alone
	  int d = h[j-1] + satBaseScore (sp, s[i], r[k]) ;
	  int v = d > f[j] ? d : f[j] ; if (e > v) v = e ;
	  hNew[j] = v ; if (v > rowMax) rowMax = v ;
	}
      int *t = h ; h = hNew ; hNew = t ;
      result = h[L] ;
      if (rowMax + sp->match * (n - 1 - i) <= threshold) { result = rowMax + sp->match * (n - 1 - i) ; break ; }
    }
  newFree (h, L+1, int) ; newFree (hNew, L+1, int) ; newFree (f, L+1, int) ;
  return result ;
}

long long satCyclicAligned = 0, satCyclicUndecided = 0, satCyclicResolved = 0 ;
long long satCyclicEarly = 0, satCyclicUpper = 0, satCyclicLower = 0 ;

#define NTRY 64 // the most last-row cells whose starts are tried when the bounds leave a pair undecided

int satCyclicScore (const char *a, int aLen, const char *b, int bLen, const SatScheme *sp, int threshold)
{
  const char *r = aLen >= bLen ? a : b, *s = aLen >= bLen ? b : a ; // r the longer, on the columns
  int L = aLen >= bLen ? aLen : bLen, n = aLen >= bLen ? bLen : aLen ;
  static const I64 MINUS_INF = -((I64)1 << 60) ;
  I64 GO = (I64)sp->gapOpen << 32, GE = (I64)sp->gapExtend << 32 ;
  I64 WRAP = L ; // added to the start label on crossing from column L-1 to 0, so labels count cycles
  I64 *h = new (L, I64), *hNew = new (L, I64), *f = new (L, I64), *gIn = new (L, I64) ;
  for (int j = 0 ; j < L ; ++j) { h[j] = (j+1 == L) ? 0 : j+1 ; f[j] = MINUS_INF ; } // score 0, start next column
  __atomic_fetch_add (&satCyclicAligned, 1, __ATOMIC_RELAXED) ;
  for (int i = 0 ; i < n ; ++i)
    { int c = s[i] ;
      for (int j = 0 ; j < L ; ++j)
	{ I64 o = h[j] + GO, x = f[j] + GE ; f[j] = o > x ? o : x ;			// s[i] alone
	  I64 d = (j ? h[j-1] : h[L-1] + WRAP) + ((I64)satBaseScore (sp, c, r[j]) << 32) ;
	  hNew[j] = d > f[j] ? d : f[j] ;
	}
      // r[j] alone, along the row and round the circle: a first pass from column 0 entering from column
      // L-1, then on round while the gap beats the one that arrived at each cell on the first pass
      I64 g = hNew[L-1] + GO + WRAP ;
      for (int j = 0 ; j < L ; ++j)
	{ gIn[j] = g ;
	  if (g > hNew[j]) { hNew[j] = g ; g += GE ; }
	  else { I64 o = hNew[j] + GO, x = g + GE ; g = o > x ? o : x ; }
	}
      g += WRAP ; // round again
      for (int j = 0 ; j < L && g > gIn[j] + WRAP ; ++j, g += GE) if (g > hNew[j]) hNew[j] = g ;
      I64 rowMax = hNew[0] ;
      for (int j = 1 ; j < L ; ++j) if (hNew[j] > rowMax) rowMax = hNew[j] ;
      I64 *t = h ; h = hNew ; hNew = t ;
      int bound = (int)(rowMax >> 32) + sp->match * (n - 1 - i) ;
      if (bound <= threshold) // even the best can no longer end above threshold
	{ __atomic_fetch_add (&satCyclicEarly, 1, __ATOMIC_RELAXED) ;
	  newFree (h, L, I64) ; newFree (hNew, L, I64) ; newFree (f, L, I64) ; newFree (gIn, L, I64) ;
	  return bound ;
	}
    }
  int ub = -(1 << 28), lb = -(1 << 28) ;
  for (int j = 0 ; j < L ; ++j)
    { int sc = (int)(h[j] >> 32), label = h[j] & 0xffffffff, start = label % L, wraps = label / L ;
      if (sc > ub) ub = sc ;
      if (start == ((j+1 == L) ? 0 : j+1) && wraps == 1 && sc > lb) lb = sc ; // exactly once round r
    }
  int result = lb ;
  if (ub <= threshold) __atomic_fetch_add (&satCyclicUpper, 1, __ATOMIC_RELAXED) ;	// the best anywhere is not above
  else if (lb > threshold) __atomic_fetch_add (&satCyclicLower, 1, __ATOMIC_RELAXED) ;	// a real cycle is above
  if (ub > threshold && lb <= threshold) // undecided: align from the starts suggested by the best cells
    { __atomic_fetch_add (&satCyclicUndecided, 1, __ATOMIC_RELAXED) ;
      // paths from neighbouring starts merge within a few rows, so a cell often holds a neighbour's start:
      // for each of the best cells, at column j, try the start j+1 that would complete a cycle there, and
      // the start it holds, until the pair is above threshold or the cells are not
      int tried[2*NTRY], nTried = 0 ;
      for (int k = 0 ; k < NTRY && result <= threshold ; ++k)
	{ int bj = -1 ;
	  for (int j = 0 ; j < L ; ++j) if (h[j] > MINUS_INF && (bj < 0 || h[j] > h[bj])) bj = j ;
	  if (bj < 0 || (int)(h[bj] >> 32) <= threshold) break ;
	  int starts[2] = { (bj+1) % L, (h[bj] & 0xffffffff) % L } ;
	  h[bj] = MINUS_INF ; // used
	  for (int q = 0 ; q < 2 && result <= threshold ; ++q)
	    { int x ; for (x = 0 ; x < nTried ; ++x) if (tried[x] == starts[q]) break ;
	      if (x < nTried) continue ;
	      tried[nTried++] = starts[q] ;
	      int sc = rotatedScore (s, n, r, L, starts[q], sp, threshold) ;
	      if (sc > result) result = sc ;
	    }
	}
      if (result > threshold) __atomic_fetch_add (&satCyclicResolved, 1, __ATOMIC_RELAXED) ;
    }
  newFree (h, L, I64) ; newFree (hNew, L, I64) ; newFree (f, L, I64) ; newFree (gIn, L, I64) ;
  return result ;
}

/******************* end of file **************/
