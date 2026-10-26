/*  File: satscore.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: scoring shared by satscan and satmatch: the weak, medium and strong schemes, the
 *   composition test for which is safe, and a global alignment of a sequence against a circular unit
 * Exported functions: see below
 * HISTORY:
 * Last edited: Oct 26 17:43 2026 (claude)
 * Created: Mon Oct 26 12:00:00 2026 (claude)
 *-------------------------------------------------------------------
 */

#ifndef SATSCORE_DEFINED
#define SATSCORE_DEFINED

// Three schemes of increasing strictness share the match score, so scores are comparable. A scheme is
// safe for a composition if the expected ungapped score per aligned base of random sequences of that
// composition is at most DANGER_T; simulations showed the gapped drift through random sequence then stays
// below about -0.20 per base in match units for all three schemes (see satscan-scoring.md).
#define MAX_SCHEME 3	// weak, medium and strong
#define DANGER_T (-1)

typedef struct {
  int match ;
  int transition ;      // A<->G, C<->T
  int transversion ;    // A<->C, A<->T, G<->C, G<->T
  int gapOpen ;         // score of the first base of a gap
  int gapExtend ;       // score of each further base of a gap
} SatScheme ;

extern const SatScheme satScheme[MAX_SCHEME] ;
extern const char satSchemeChar[MAX_SCHEME] ; // 'w', 'm', 's'

static inline int satBaseScore (const SatScheme *sp, int a, int b) // a, b in 0..3, ACGT
{ return a == b ? sp->match : ((a ^ b) == 2 ? sp->transition : sp->transversion) ; }

// the expected score per base of a base from composition p against one from composition q, both in ACGT
// order and summing to 1
double satExpectedScore (const SatScheme *sp, const double *p, const double *q) ;

// the least strict scheme that is safe for a sequence against random sequence of its own composition,
// or MAX_SCHEME if none is; seq is in 0..3
int satSchemeLevel (const char *seq, int len) ;

// The best score of a global alignment of all of t against the circular unit q, which may start at any
// position of q and wrap round any number of times, under scheme sp with affine gaps; t and q are in 0..3.
// The alignment stops as soon as it cannot end above threshold, as a row's best score can rise by at most
// sp->match per remaining base of t, and then returns a bound that is at most threshold.
int satGlobalScore (const char *t, int tLen, const char *q, int qLen, const SatScheme *sp, int threshold) ;

#endif

/******************* end of file **************/
