/*  File: satscan.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: statellite scanning package
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 20 23:59 2026 (claude)
 * Created: Sun Aug 16 14:02:47 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "onerep.h"
#include "array.h"  // contains utils.h
#include "seqio.h"
#include "ONElib.h"
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>

// DP values pack the alignment score into the upper 32 bits (signed) and the number of exact matches
// into the lower 32 bits. Penalties are multiples of ONE so they leave the match count alone, and
// comparing packed values compares scores first and then prefers more exact matches.
// NB the match count would carry into the score at 2^32 exact matches in one alignment.
#define ONE ((I64)1 << 32)

// The copy parse (approach A'): every DP cell carries the index of the CopyNode for the start of
// the copy of the unit that its path is in. A new node is made each time a path wraps round from the
// end of the unit to the start, so following prev from the best end recovers the copy boundaries.
typedef struct {
  U32  pos ;	   // start of a copy, relative to match->start
  U32  prev ;	   // node of the previous copy start; node 0 is the match start
} CopyNode ;

typedef struct MatchStruct {
  U64  start ;	   // position in current sequence
  I64  bestScore ; // best score so far at end of a query unit
  I64  maxScore ;  // best cell score so far, for the x-drop
  U64  bestEnd ;   // end position of bestScore
  I64 *score ;     // length q->len - must be positive
  I64 *gap ;       // length q->len
  U32 *copy ;      // length q->len: CopyNode index for each score cell
  U32 *gapCopy ;   // length q->len: CopyNode index for each gap cell
  Array aNode ;    // of CopyNode
  U32  nodeLimit ; // compact aNode when it grows past this
  U32  bestNode ;  // CopyNode of the copy that bestScore ends in
  U32 *parse ;     // internal copy boundaries relative to start, made by matchFinish() if reported
  int  nParse ;
  int  query ;     // index in aQuery
  int  target ;    // index of the target sequence, in file order
  U8   schemes ;   // bit s set if scheme s was used somewhere in [start, bestEnd), set if reported
  struct MatchStruct  *next ;      // index of nextMatch in aMatch
} Match ;

typedef struct {
  int match ;
  int transition ;      /* A<->G, C<->T */
  int transversion ;    /* A<->C, A<->T, G<->C, G<->T */
  I64 gapOpen ;         /* gap penalty per position, times ONE */
  I64 gapExtend ;       /* gap extend per position, times ONE */
} ScoreParams ;

// By default three scoring schemes of increasing strictness are used. For each query strand and each
// region of the target, satscan uses the least strict scheme under which chance matches are unlikely
// there, judged from base composition (see dangerLevels()). A scheme given with -score is used
// everywhere instead. The schemes share the match score, so scores, -smin and the x-drop are comparable.
#define MAX_SCHEME 3
static ScoreParams defaultScheme[MAX_SCHEME] = {
  { 2, -1, -3, -15*ONE, -1*ONE },	// weak
  { 2, -2, -4, -20*ONE, -2*ONE },	// medium
  { 2, -3, -6, -25*ONE, -2*ONE } } ;	// strong
static char schemeChar[MAX_SCHEME] = { 'w', 'm', 's' } ; // for the output; a -score scheme is 'u'

// A region is dangerous for a query strand and scheme if the expected ungapped score per aligned base
// of the query against random sequence of the region's composition is above DANGER_T (score units).
// Where it is at or below DANGER_T, simulations showed the gapped drift through random sequence stays
// below about -0.20 per base in match units for all three schemes (see satscan-scoring.md).
#define DANGER_T (-1)

typedef struct {
  int   len ;
  int   isRC ;
  int   index ;     // 1-based position of the entry in the query file, counting any that were skipped
  char *name ;      // shared between a forward query and its reverse complement
  char *seq ;
  I64   startScore[MAX_SCHEME][32][4] ;
  I64   endScore[MAX_SCHEME][32][4] ;
  I64   xDrop ;     // cells more than this below the match's best cell so far are dropped, times ONE
  I64   dangerWeight[MAX_SCHEME][4] ; // len times (expected score of target base a vs a random query base - DANGER_T)
  U64   dangerMinLen ;		// dangerous regions shorter than this are ignored
} Query ;

typedef struct {		// a region of the target scanned with scheme level, or walled if level == nScheme
  U64  start, end ;
  int  level ;
} LevelRun ;

typedef struct { U64 start, end ; } Segment ;
typedef struct { U64 pos ; int delta ; } Event ;	// a dangerous region starts (+1) or ends (-1)

// We are going to require an ungapped positive 32bp match to the start of the query
// to start, and for symmetry an ungapped positive 32bp match to the end of the query
// to end. Extend by affine DP in between.

void queryDestroy (Array a)
{
  for (int i = 0 ; i < arrayMax(a) ; ++i)
    { Query *q = arrp(a,i,Query) ;
      newFree (q->seq, q->len, char) ;
      if (!q->isRC) newFree (q->name, strlen(q->name)+1, char) ;
    }
  arrayDestroy (a) ;
}

static I64 baseScore[MAX_SCHEME][4][4] ;	// score times ONE, plus 1 for an exact match

static void makeQueryScores (Query *q, int nScheme)
{
  for (int s = 0 ; s < nScheme ; ++s)
    for (int i = 0 ; i < 32 ; ++i)
      for (int j = 0 ; j < 4 ; ++j)
	{ q->startScore[s][i][j] = baseScore[s][(int)q->seq[i%q->len]][j] ;
	  q->endScore[s][i][j] = baseScore[s][(int)q->seq[(32*q->len+i-32)%q->len]][j] ;
	}
  I64 n[4] = { 0, 0, 0, 0 } ;
  for (int k = 0 ; k < q->len ; ++k) ++n[(int)q->seq[k]] ;
  for (int s = 0 ; s < nScheme ; ++s)
    for (int a = 0 ; a < 4 ; ++a)
      { I64 w = -(I64)DANGER_T * q->len ;
	for (int b = 0 ; b < 4 ; ++b) w += n[b] * (baseScore[s][a][b] >> 32) ;
	q->dangerWeight[s][a] = w ;
      }
  q->dangerMinLen = q->len > 64 ? q->len : 64 ;
}

typedef struct {		// read-only once scanning starts, so can be shared between threads
  Array       aQuery ;
  ScoreParams sp[MAX_SCHEME] ;
  int         nScheme ;		// 1 if -score was given, else MAX_SCHEME
  U64         scoreMin ;
} ScanParams ;

typedef struct {		// everything that is private to one scanning thread
  Match **active ;		// [nQuery] heads of the lists of active matches for each query
  Array   aMatch ;		// Match* to report: score and gap are freed, so these are small
  I64    *gapIn ;		// [maxLen] scratch: right gap value arriving at each cell on the first pass
  int     maxLen ;		// longest query
  U64     totStart, totMatch, totMatchLen, totSameEnd ; // summed over all the targets this worker did
  Array  *aRun ;		// [nQuery] LevelRun for the current target, sorted, level > 0 only
  U64    *runCursor ;		// [nQuery] first run that can still matter as the scan moves along
  Array   aSeg[MAX_SCHEME] ;	// scratch Segment lists for dangerLevels()
  Array   aEvent ;		// scratch for dangerLevels()
  U64    *regionN, *regionLen ; // [nQuery*(MAX_SCHEME+1)] count and length of regions at each level
} Worker ;

typedef struct {
  int   target ;		// index of this target sequence, in file order
  char *name ;
  char *seq ;			// owned by the job if and only if ownsSeq
  U64   len ;
  bool  ownsSeq ;		// true if seq is a copy that the worker must newFree
} Job ;

typedef struct {
  U64    nStart, nMatch, matchLen ; // nMatch and matchLen are after removeSameEnd()
  U64    nSameEnd ;		// matches removed because a better one of the same query ends at the same place
  double seconds ;		// wall time
} ScanStats ;

static double wallTime (void)
{ struct timespec t ;
  clock_gettime (CLOCK_MONOTONIC, &t) ;
  return t.tv_sec + 1e-9*t.tv_nsec ;
}

static Worker *workerCreate (ScanParams *scan)
{ Worker *w = new0 (1, Worker) ;
  w->active = new0 (arrayMax(scan->aQuery), Match*) ;
  w->aMatch = arrayCreate (256, Match*) ;
  for (int j = 0 ; j < arrayMax(scan->aQuery) ; ++j)
    if (arrp(scan->aQuery, j, Query)->len > w->maxLen) w->maxLen = arrp(scan->aQuery, j, Query)->len ;
  w->gapIn = new (w->maxLen, I64) ;
  int nQuery = arrayMax(scan->aQuery) ;
  w->aRun = new (nQuery, Array) ;
  for (int j = 0 ; j < nQuery ; ++j) w->aRun[j] = arrayCreate (64, LevelRun) ;
  w->runCursor = new0 (nQuery, U64) ;
  for (int s = 0 ; s < MAX_SCHEME ; ++s) w->aSeg[s] = arrayCreate (1024, Segment) ;
  w->aEvent = arrayCreate (1024, Event) ;
  w->regionN = new0 (nQuery*(MAX_SCHEME+1), U64) ;
  w->regionLen = new0 (nQuery*(MAX_SCHEME+1), U64) ;
  return w ;
}

static void workerDestroy (Worker *w, ScanParams *scan)
{ int nQuery = arrayMax(scan->aQuery) ;
  newFree (w->active, nQuery, Match*) ;
  newFree (w->gapIn, w->maxLen, I64) ;
  for (int j = 0 ; j < nQuery ; ++j) arrayDestroy (w->aRun[j]) ;
  newFree (w->aRun, nQuery, Array) ;
  newFree (w->runCursor, nQuery, U64) ;
  for (int s = 0 ; s < MAX_SCHEME ; ++s) arrayDestroy (w->aSeg[s]) ;
  arrayDestroy (w->aEvent) ;
  newFree (w->regionN, nQuery*(MAX_SCHEME+1), U64) ;
  newFree (w->regionLen, nQuery*(MAX_SCHEME+1), U64) ;
  arrayDestroy (w->aMatch) ;
  newFree (w, 1, Worker) ;
}

static inline U32 copyNodeAdd (Match *match, U64 pos, U32 prev) // pos is in target coordinates
{
  U32 n = arrayMax(match->aNode) ;
  CopyNode *c = arrayp(match->aNode, n, CopyNode) ;
  c->pos = pos - match->start ; c->prev = prev ;
  return n ;
}

// Keep only the nodes reachable from live cells and bestNode, renumbering them in order. A node's prev
// is always made before it, so renumbering in order keeps prev < node, and node 0 stays the start.
static void copyNodeCompact (Match *match, int len)
{
  U32  n = arrayMax(match->aNode) ;
  CopyNode *node = arrp(match->aNode, 0, CopyNode) ;
  U32 *map = new0 (n, U32) ; // first marks the reachable nodes with 1, then holds the new index
  map[0] = 1 ;
  for (int k = 0 ; k < len ; ++k)
    { if (match->score[k] < ONE) match->copy[k] = 0 ;
      if (match->gap[k] < ONE) match->gapCopy[k] = 0 ;
      for (U32 c = match->copy[k] ; !map[c] ; c = node[c].prev) map[c] = 1 ;
      for (U32 c = match->gapCopy[k] ; !map[c] ; c = node[c].prev) map[c] = 1 ;
    }
  for (U32 c = match->bestNode ; !map[c] ; c = node[c].prev) map[c] = 1 ;
  U32 nNew = 0 ;
  for (U32 c = 0 ; c < n ; ++c)
    if (map[c])
      { node[nNew].pos = node[c].pos ; node[nNew].prev = map[node[c].prev] ; // prev < c so already mapped
	map[c] = nNew++ ;
      }
  arrayMax(match->aNode) = nNew ;
  for (int k = 0 ; k < len ; ++k)
    { match->copy[k] = map[match->copy[k]] ; match->gapCopy[k] = map[match->gapCopy[k]] ; }
  match->bestNode = map[match->bestNode] ;
  newFree (map, n, U32) ;
  match->nodeLimit = 2*nNew + 4*len ;
}

// fill match->parse with the internal copy boundaries of the best alignment, in increasing order
static void matchParse (Match *match, Query *q)
{
  CopyNode *node = arrp(match->aNode, 0, CopyNode) ;
  int n = 0 ;
  for (U32 c = match->bestNode ; c ; c = node[c].prev) ++n ;
  U64 endStart = match->bestEnd - 32 - match->start ; // the final 32bp, aligned ungapped to the end of the unit
  for (int k = 0 ; k < 32 ; ++k) if ((32 - k) % q->len == 0) ++n ; // copies starting inside that
  match->nParse = n ;
  match->parse = new (n ? n : 1, U32) ;
  int i = n ;
  for (int k = 31 ; k >= 0 ; --k) if ((32 - k) % q->len == 0) match->parse[--i] = endStart + k ;
  for (U32 c = match->bestNode ; c ; c = node[c].prev) match->parse[--i] = node[c].pos ;
}

/********** dangerous regions, where chance matches are likely, and the scheme to use in them **********/

// positive root of sum_a g_a exp(lambda W_a) = 1, which exists when sum_a g_a W_a < 0 and some W_a > 0
static double dangerLambda (I64 *W, double *g)
{
  double lo = 0, hi = 1e-6 ;
  for (;;)
    { double f = -1 ; for (int a = 0 ; a < 4 ; ++a) f += g[a] * exp (hi * W[a]) ;
      if (f > 0) break ;
      hi *= 2 ;
    }
  for (int it = 0 ; it < 100 ; ++it)
    { double mid = 0.5 * (lo + hi), f = -1 ;
      for (int a = 0 ; a < 4 ; ++a) f += g[a] * exp (mid * W[a]) ;
      if (f > 0) hi = mid ; else lo = mid ;
    }
  return 0.5 * (lo + hi) ;
}

// Add to aSeg the segments of seq[from,to) whose summed weight is at least sStar and that are at least
// minLen long, by the one-pass reset-at-zero method.
static void dangerSegments (char *seq, U64 from, U64 to, I64 *W, double sStar, U64 minLen, Array aSeg)
{
  I64 sum = 0, maxSum = 0 ;
  U64 start = from, maxEnd = from ;
  for (U64 p = from ; p < to ; ++p)
    { int a = seq[p] ;
      if ((sum += W[a]) > 0)
	{ if (sum > maxSum) { maxSum = sum ; maxEnd = p+1 ; } continue ; }
      if (maxSum >= sStar && maxEnd - start >= minLen)
	{ Segment *sg = arrayp (aSeg, arrayMax(aSeg), Segment) ; sg->start = start ; sg->end = maxEnd ; }
      sum = maxSum = 0 ; start = p+1 ;
    }
  if (maxSum >= sStar && maxEnd - start >= minLen)
    { Segment *sg = arrayp (aSeg, arrayMax(aSeg), Segment) ; sg->start = start ; sg->end = maxEnd ; }
}

static int eventOrder (const void *a, const void *b)
{ const Event *x = (const Event*)a, *y = (const Event*)b ;
  return x->pos < y->pos ? -1 : x->pos > y->pos ? 1 : 0 ;
}

// For each query, fill w->aRun[j] with the regions of the target where the weak scheme, then each
// stricter one in turn, is dangerous. A region dangerous for schemes 0..s-1 but not s is scanned with
// scheme s; one dangerous for all of them is walled: no match may start in it or continue through it.
// Regions for scheme s are looked for only inside those for scheme s-1, as the schemes get stricter.
// seq is in 0..3 (dna2index4Conv, so N counts as A); a scheme given by -score is never masked.
static void dangerLevels (ScanParams *scan, Worker *w, char *seq, U64 len)
{
  int nQuery = arrayMax(scan->aQuery) ;
  for (int j = 0 ; j < nQuery ; ++j) { arrayMax(w->aRun[j]) = 0 ; w->runCursor[j] = 0 ; }
  if (scan->nScheme == 1 || !len) return ;
  double g[4] = { 0, 0, 0, 0 } ;
  for (U64 p = 0 ; p < len ; ++p) ++g[(int)seq[p]] ;
  for (int a = 0 ; a < 4 ; ++a) g[a] /= len ;
  for (int j = 0 ; j < nQuery ; ++j)
    { Query *q = arrp(scan->aQuery, j, Query) ;
      Array prev = 0 ; // 0 means the whole target
      for (int s = 0 ; s < scan->nScheme ; ++s)
	{ Array aSeg = w->aSeg[s] ; arrayMax(aSeg) = 0 ;
	  if (prev && !arrayMax(prev)) { prev = aSeg ; continue ; } // nothing dangerous for a weaker scheme
	  I64 *W = q->dangerWeight[s], wMax = W[0] ;
	  double mean = 0 ;
	  for (int a = 0 ; a < 4 ; ++a) { mean += g[a] * W[a] ; if (W[a] > wMax) wMax = W[a] ; }
	  if (mean >= 0) // the composition of the whole target is dangerous: all of prev is
	    { if (!prev) { Segment *sg = arrayp (aSeg, 0, Segment) ; sg->start = 0 ; sg->end = len ; }
	      else for (U64 k = 0 ; k < arrayMax(prev) ; ++k) array (aSeg, k, Segment) = arr (prev, k, Segment) ;
	    }
	  else if (wMax > 0)
	    { double sStar = log (1e8) / dangerLambda (W, g) ; // about one chance region per 100Mb
	      if (!prev) dangerSegments (seq, 0, len, W, sStar, q->dangerMinLen, aSeg) ;
	      else
		for (U64 k = 0 ; k < arrayMax(prev) ; ++k)
		  dangerSegments (seq, arrp(prev,k,Segment)->start, arrp(prev,k,Segment)->end,
				  W, sStar, q->dangerMinLen, aSeg) ;
	    }
	  prev = aSeg ;
	}
      // nested segments to runs: the level at a position is the number of schemes it is dangerous for
      Array aEvent = w->aEvent ; arrayMax(aEvent) = 0 ;
      for (int s = 0 ; s < scan->nScheme ; ++s)
	for (U64 k = 0 ; k < arrayMax(w->aSeg[s]) ; ++k)
	  { Segment *sg = arrp(w->aSeg[s], k, Segment) ;
	    Event *e = arrayp (aEvent, arrayMax(aEvent), Event) ; e->pos = sg->start ; e->delta = 1 ;
	    e = arrayp (aEvent, arrayMax(aEvent), Event) ; e->pos = sg->end ; e->delta = -1 ;
	  }
      arraySort (aEvent, eventOrder) ;
      int level = 0 ;
      for (U64 k = 0 ; k < arrayMax(aEvent) ; )
	{ U64 pos = arrp(aEvent,k,Event)->pos ;
	  while (k < arrayMax(aEvent) && arrp(aEvent,k,Event)->pos == pos) level += arrp(aEvent,k++,Event)->delta ;
	  if (k < arrayMax(aEvent) && level > 0)
	    { U64 end = arrp(aEvent,k,Event)->pos ;
	      Array aRun = w->aRun[j] ;
	      LevelRun *last = arrayMax(aRun) ? arrp(aRun, arrayMax(aRun)-1, LevelRun) : 0 ;
	      if (last && last->end == pos && last->level == level) last->end = end ;
	      else { LevelRun *r = arrayp (aRun, arrayMax(aRun), LevelRun) ; r->start = pos ; r->end = end ; r->level = level ; }
	    }
	}
      for (U64 k = 0 ; k < arrayMax(w->aRun[j]) ; ++k)
	{ LevelRun *r = arrp(w->aRun[j], k, LevelRun) ;
	  ++w->regionN[j*(MAX_SCHEME+1) + r->level] ; w->regionLen[j*(MAX_SCHEME+1) + r->level] += r->end - r->start ;
	}
    }
}

// the levels of the 65 positions p-32..p+32 for query j, advancing its cursor; returns true if they
// are all the same, in which case only lv[0] is set
static inline bool windowLevels (Worker *w, int j, U64 p, int *lv)
{
  Array aRun = w->aRun[j] ;
  U64 nRun = arrayMax(aRun), r = w->runCursor[j] ;
  LevelRun *run = nRun ? arrp(aRun, 0, LevelRun) : 0 ;
  while (r < nRun && run[r].end <= p - 32) ++r ;
  w->runCursor[j] = r ;
  if (r == nRun || run[r].start > p + 32) { lv[0] = 0 ; return true ; }
  if (run[r].start <= p - 32 && run[r].end > p + 32) { lv[0] = run[r].level ; return true ; }
  for (int k = 0 ; k < 65 ; ++k)
    { U64 x = p - 32 + k ;
      while (r < nRun && run[r].end <= x) ++r ;
      lv[k] = (r < nRun && run[r].start <= x) ? run[r].level : 0 ;
    }
  return false ;
}

// the bit set of levels used over [start, end) for query j
static U8 schemesUsed (Worker *w, int j, U64 start, U64 end)
{
  Array aRun = w->aRun[j] ;
  U8 bits = 0 ; U64 covered = 0 ;
  for (U64 k = 0 ; k < arrayMax(aRun) ; ++k)
    { LevelRun *r = arrp(aRun, k, LevelRun) ;
      if (r->end <= start) continue ;
      if (r->start >= end) break ;
      U64 lo = r->start > start ? r->start : start, hi = r->end < end ? r->end : end ;
      bits |= 1 << r->level ; covered += hi - lo ;
    }
  if (covered < end - start) bits |= 1 ;
  return bits ;
}

// match is done and unlinked: keep it for reporting if it scores well enough, else free it
static void matchFinish (Worker *w, ScanParams *scan, Match *match, Query *q, ScanStats *st)
{
  newFree (match->score, q->len, I64) ; match->score = 0 ;
  newFree (match->gap, q->len, I64) ; match->gap = 0 ;
  newFree (match->copy, q->len, U32) ; match->copy = 0 ;
  newFree (match->gapCopy, q->len, U32) ; match->gapCopy = 0 ;
  bool isReport = (match->bestScore >> 32) > (I64)scan->scoreMin ;
  if (isReport) matchParse (match, q) ;
  arrayDestroy (match->aNode) ; match->aNode = 0 ;
  if (isReport)
    { match->next = 0 ;
      match->schemes = schemesUsed (w, match->query, match->start, match->bestEnd) ;
      array(w->aMatch, arrayMax(w->aMatch), Match*) = match ;
      ++st->nMatch ; st->matchLen += match->bestEnd - match->start ;
    }
  else
    newFree (match, 1, Match) ;
}

// groups matches with the same target, query (so strand) and end, best first: higher score, then
// higher match count, then earlier start, so the choice is deterministic
static int sameEndOrder (const void *a, const void *b)
{
  const Match *x = *(Match* const*)a, *y = *(Match* const*)b ;
  if (x->target != y->target) return x->target < y->target ? -1 : 1 ;
  if (x->query != y->query) return x->query < y->query ? -1 : 1 ;
  if (x->bestEnd != y->bestEnd) return x->bestEnd < y->bestEnd ? -1 : 1 ;
  if (x->bestScore != y->bestScore) return x->bestScore > y->bestScore ? -1 : 1 ; // packed: score, then count
  if (x->start != y->start) return x->start < y->start ? -1 : 1 ;
  return 0 ;
}

// keep only the best of the matches in aMatch from i0 on (all from one target) that share a query and
// end, updating the counts in st
static void removeSameEnd (Array aMatch, U64 i0, ScanStats *st)
{
  U64 nIn = arrayMax(aMatch) - i0 ;
  if (nIn < 2) return ;
  Match **mm = arrp(aMatch, i0, Match*) ;
  qsort (mm, nIn, sizeof(Match*), sameEndOrder) ;
  U64 n = 0 ;
  for (U64 i = 0 ; i < nIn ; ++i)
    { Match *m = mm[i], *keep = n ? mm[n-1] : 0 ;
      if (keep && keep->query == m->query && keep->bestEnd == m->bestEnd)
	{ --st->nMatch ; st->matchLen -= m->bestEnd - m->start ; ++st->nSameEnd ;
	  newFree (m->parse, m->nParse ? m->nParse : 1, U32) ; newFree (m, 1, Match) ;
	}
      else
	mm[n++] = m ;
    }
  arrayMax(aMatch) = i0 + n ;
}


static void scanTarget (ScanParams *scan, Worker *w, Job *job, ScanStats *st)
{
  Array aQuery = scan->aQuery ;
  int nQuery = arrayMax(aQuery), nScheme = scan->nScheme ;
  memset (st, 0, sizeof(ScanStats)) ;
  U64 i0 = arrayMax(w->aMatch) ; // this target's matches are added from here
  double t0 = wallTime () ;

  dangerLevels (scan, w, job->seq, job->len) ; // find where each query needs a stricter scheme

  char *s = job->seq + 32 ;
  for (U64 i = 32 ; i < job->len ; ++i, ++s)
    { for (int j = 0 ; j < nQuery ; ++j)
	{ Query *q  = arrp(aQuery, j, Query) ;
	  char  *qs = q->seq ;
	  int lv[65] ; // scheme levels at target positions i-32..i+32, so lv[32] is at i, unless uniform
	  bool isUniform = windowLevels (w, j, i, lv) ;
	  int li = isUniform ? lv[0] : lv[32] ;
	  if (li == nScheme) // a wall: end the matches here, and start none
	    { Match *match = w->active[j], *next ;
	      for ( ; match ; match = next) { next = match->next ; matchFinish (w, scan, match, q, st) ; }
	      w->active[j] = 0 ;
	      continue ;
	    }
	  const ScoreParams *sp = &scan->sp[li] ;
	  I64 *bs = baseScore[li][(int)*s] ;
	  // a start needs t[i-32..i-1] to score positively ungapped against q[0..31], and an end needs
	  // t[i+1..i+32] against the end of the unit; each base is scored with the scheme at its position
	  I64 newStart = 0, newEnd = 0 ; // the end window starts at t[i+1], after the cell for q[len-33] at t[i]
	  bool isEndOK = i < job->len - 32 ;
	  if (isUniform)
	    { for (int k = 0 ; k < 32 ; ++k) newStart += q->startScore[li][k][(int)s[k-32]] ;
	      if (isEndOK) for (int k = 0 ; k < 32 ; ++k) newEnd += q->endScore[li][k][(int)s[k+1]] ;
	    }
	  else
	    { bool isStartOK = true ;
	      for (int k = 0 ; k < 32 ; ++k) if (lv[k] == nScheme) isStartOK = false ;
	      if (isStartOK) for (int k = 0 ; k < 32 ; ++k) newStart += q->startScore[lv[k]][k][(int)s[k-32]] ;
	      for (int k = 33 ; k < 65 ; ++k) if (lv[k] == nScheme) isEndOK = false ;
	      if (isEndOK) for (int k = 0 ; k < 32 ; ++k) newEnd += q->endScore[lv[33+k]][k][(int)s[k+1]] ;
	    }
	  if (newStart >= ONE) newStart += bs[(int)qs[32 % q->len]] ; // the seed cell is (q[32%len], t[i])
	  Match *lastMatch = 0, *next ;
	  for (Match *match = w->active[j] ; match ; match = next)
	    {
	      next = match->next ; // save before match can be freed below
	      // first rotate: the cell moving from the end of the unit to the start begins a new copy at i
	      I64 *ms = match->score ;
	      I64 msn = ms[q->len-1] ; memmove (ms+1, ms, (q->len-1)*sizeof(I64)) ; ms[0] = msn ;
	      U32 *mc = match->copy ;
	      U32 mcn = mc[q->len-1] ; memmove (mc+1, mc, (q->len-1)*sizeof(U32)) ;
	      mc[0] = msn ? copyNodeAdd (match, i, mcn) : 0 ;
	      // then add diagonal match scores
	      for (int k = 0 ; k < q->len ; ++k, ++ms) if (*ms) *ms += bs[(int)qs[k]] ;
	      // next right gaps. These delete query bases, so when a gap wraps from the end of the unit
	      // to the start, the next copy starts at i+1. Only make its CopyNode if the gap is used.
	      ms = match->score ; mc = match->copy ;
	      I64 d = ms[q->len-1] + sp->gapOpen ;
	      U32 dc = mc[q->len-1] ;
	      bool isWrap = true ; // d has wrapped since dc was set
	      I64 *gi = w->gapIn ;
	      for (int k = 0 ; k < q->len ; ++k, ++ms, ++mc)
		if ((*gi++ = d) > *ms)
		  { if (isWrap) { dc = copyNodeAdd (match, i+1, dc) ; isWrap = false ; }
		    *ms = d ; *mc = dc ; d += sp->gapExtend ;
		  }
		else // extend, or open from here
		  { I64 e = d + sp->gapExtend, o = *ms + sp->gapOpen ;
		    if (o > e) { d = o ; dc = *mc ; isWrap = false ; } else d = e ;
		  }
	      // Extend the wrap-around while it beats the gap that arrived at each cell on the first pass:
	      // it can pass cells it doesn't improve, but once below that value it can't win further on.
	      ms = match->score ; mc = match->copy ; gi = w->gapIn ;
	      bool isWrap2 = true ; // d has wrapped again
	      for (int k = 0 ; k < q->len && d > gi[k] ; ++k, d += sp->gapExtend)
		if (d > ms[k])
		  { if (isWrap2)
		      { if (isWrap) dc = copyNodeAdd (match, i+1, dc) ; // gap went right round the unit unused
			dc = copyNodeAdd (match, i+1, dc) ; isWrap2 = false ;
		      }
		    ms[k] = d ; mc[k] = dc ;
		  }
	      // next the down gaps, always: a gap can open from any cell, even after a low patch
	      { ms = match->score ; mc = match->copy ;
		I64 *mg = match->gap ;
		U32 *mgc = match->gapCopy ;
		for (int k = 0 ; k < q->len ; ++k, ++ms, ++mg, ++mc, ++mgc)
		  if (*mg > *ms) { *ms = *mg ; *mc = *mgc ; *mg += sp->gapExtend ; }
		  else // extend, or open from here
		    { I64 e = *mg + sp->gapExtend, o = *ms + sp->gapOpen ;
		      if (o > e) { *mg = o ; *mgc = *mc ; } else *mg = e ;
		    }
	      }
	      if (arrayMax(match->aNode) > match->nodeLimit) copyNodeCompact (match, q->len) ;
	      // finally, check if done, and if so then process match
	      ms = match->score ;
	      // A cell is alive if its score is > 0, whatever the match count, and it is within xDrop of
	      // the best cell so far, so no reported path can lose more than xDrop over any stretch.
	      bool isDone = true ;
	      I64 minLive = match->maxScore - q->xDrop, colMax = 0 ;
	      if (minLive < ONE) minLive = ONE ;
	      for (int k = 0 ; k < q->len ; ++k, ++ms)
		if (*ms >= minLive) { isDone = false ; if (*ms > colMax) colMax = *ms ; }
		else *ms = 0 ;
	      if (colMax > match->maxScore) match->maxScore = colMax ;
	      if (isDone)
		{ if (lastMatch) lastMatch->next = next ; else w->active[j] = next ;
		  matchFinish (w, scan, match, q, st) ;
		  continue ; // match is gone: skip the code below and don't set lastMatch
		}
	      else
		{
		  // update bestScore if necessary
		  int kEnd = (33*(q->len-1)) % q->len ; // q[-33 mod len]: the end window covers the next 32
		  I64 mEnd = match->score[kEnd] ;
		  if ((newEnd >= ONE) && mEnd && (mEnd + newEnd > match->bestScore))
		    { match->bestScore = mEnd + newEnd ; match->bestEnd = i + 33 ;
		      match->bestNode = match->copy[kEnd] ;
		    }
		}

	      if ((match->score[32 % q->len] >> 32) >= (newStart >> 32)) newStart = 0 ; // don't restart a new one: compare scores only, not match counts

	      lastMatch = match ;
	    }

	  if (newStart >= ONE)
	    { ++st->nStart ;
	      Match *match = new0 (1, Match) ;
	      match->query = j ; match->target = job->target ; match->start = i-32 ;
	      match->next = w->active[j] ; w->active[j] = match ; // link into list
	      I64 *ms = match->score = new0 (q->len, I64) ;
	      match->gap = new0 (q->len, I64) ;
	      U32 *mc = match->copy = new0 (q->len, U32) ;
	      match->gapCopy = new0 (q->len, U32) ;
	      match->aNode = arrayCreate (1024, CopyNode) ;
	      match->nodeLimit = 4*q->len + 1024 ;
	      U32 c = copyNodeAdd (match, i-32, 0) ; // node 0 is the match start
	      for (int m = q->len ; m <= 32 ; m += q->len) c = copyNodeAdd (match, i-32+m, c) ; // copies in the seed
	      int k = 32 % q->len ;
	      ms[k] = newStart ; mc[k] = c ;
	      match->maxScore = newStart ;
	      I64 d = newStart + sp->gapOpen ;
	      if (d >= ONE)
		{ match->gap[k] = d ; match->gapCopy[k] = c ;
		  while (d >= ONE)
		    { if (++k == q->len) { k = 0 ; c = copyNodeAdd (match, i+1, c) ; }
		      if (ms[k] >= d) break ;
		      ms[k] = d ; mc[k] = c ;
		      d += sp->gapExtend ;
		    }
		}
	    }
	}
    }

  for (int j = 0 ; j < nQuery ; ++j) // flush matches still active at end of sequence
    { Query *q = arrp(aQuery, j, Query) ;
      Match *match = w->active[j], *next ;
      for ( ; match ; match = next)
	{ next = match->next ;
	  matchFinish (w, scan, match, q, st) ;
	}
      w->active[j] = 0 ;
    }
  removeSameEnd (w->aMatch, i0, st) ;
  st->seconds = wallTime () - t0 ;
}

static void reportTarget (Job *job, ScanStats *st)
{
  fprintf (stderr, "%s length %llu: %llu starts, %llu matches, total match length %llu, %.2f s\n",
	   job->name, job->len, st->nStart, st->nMatch, st->matchLen, st->seconds) ;
}

// scan one target, report it to stderr, update the worker's totals and release the job's sequence
static void workerRun (ScanParams *scan, Worker *w, Job *job)
{
  ScanStats st ;
  scanTarget (scan, w, job, &st) ;
  reportTarget (job, &st) ;
  w->totStart += st.nStart ; w->totMatch += st.nMatch ; w->totMatchLen += st.matchLen ; w->totSameEnd += st.nSameEnd ;
  if (job->ownsSeq) newFree (job->seq, job->len, char) ;
}

// Bounded queue of jobs from the master to the worker threads. Any idle worker takes the next job,
// so a long target never holds up others. Bounding it limits memory to about 2 sequences per thread.

typedef struct {
  pthread_mutex_t lock ;
  pthread_cond_t  notEmpty, notFull ;
  Job            *jobs ;		// ring buffer
  int             capacity, head, count ;
  bool            isClosed ;		// no more jobs will be pushed
} JobQueue ;

static void queueInit (JobQueue *q, int capacity)
{
  pthread_mutex_init (&q->lock, 0) ;
  pthread_cond_init (&q->notEmpty, 0) ; pthread_cond_init (&q->notFull, 0) ;
  q->jobs = new (capacity, Job) ;
  q->capacity = capacity ; q->head = q->count = 0 ; q->isClosed = false ;
}

static void queueDestroy (JobQueue *q)
{
  pthread_mutex_destroy (&q->lock) ;
  pthread_cond_destroy (&q->notEmpty) ; pthread_cond_destroy (&q->notFull) ;
  newFree (q->jobs, q->capacity, Job) ;
}

static void queuePush (JobQueue *q, Job *job) // blocks while the queue is full
{
  pthread_mutex_lock (&q->lock) ;
  while (q->count == q->capacity) pthread_cond_wait (&q->notFull, &q->lock) ;
  q->jobs[(q->head + q->count++) % q->capacity] = *job ;
  pthread_cond_signal (&q->notEmpty) ;
  pthread_mutex_unlock (&q->lock) ;
}

static bool queuePop (JobQueue *q, Job *job) // blocks while empty; false when closed and empty
{
  pthread_mutex_lock (&q->lock) ;
  while (!q->count && !q->isClosed) pthread_cond_wait (&q->notEmpty, &q->lock) ;
  if (!q->count) { pthread_mutex_unlock (&q->lock) ; return false ; }
  *job = q->jobs[q->head] ;
  q->head = (q->head + 1) % q->capacity ; --q->count ;
  pthread_cond_signal (&q->notFull) ;
  pthread_mutex_unlock (&q->lock) ;
  return true ;
}

static void queueClose (JobQueue *q)
{
  pthread_mutex_lock (&q->lock) ;
  q->isClosed = true ;
  pthread_cond_broadcast (&q->notEmpty) ;
  pthread_mutex_unlock (&q->lock) ;
}

typedef struct {
  ScanParams *scan ;
  JobQueue   *queue ;
  Worker     *worker ;
  pthread_t   thread ;
} ThreadArg ;

static void *workerMain (void *arg)
{
  ThreadArg *ta = (ThreadArg*) arg ;
  Job job ;
  while (queuePop (ta->queue, &job)) workerRun (ta->scan, ta->worker, &job) ;
  return 0 ;
}

static int nQueryFwd ; // number of forward queries, for matchOrder()

// sorts by target in file order, then start position, then query with each forward query before its reverse complement
static int matchOrder (const void *a, const void *b)
{
  const Match *x = *(Match* const*)a, *y = *(Match* const*)b ;
  if (x->target != y->target) return x->target < y->target ? -1 : 1 ;
  if (x->start != y->start) return x->start < y->start ? -1 : 1 ;
  int xq = x->query % nQueryFwd, yq = y->query % nQueryFwd ;
  if (xq != yq) return xq < yq ? -1 : 1 ;
  if (x->query != y->query) return x->query < y->query ? -1 : 1 ;
  return 0 ;
}

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: satscan <options>* <sat_seqfile> <genome_seqfile>\n"
    "  Finds tandem matches of complete units with affine gap alignment. By default uses a hierarchical\n"
    "  scoring scheme to avoid composition-based false positives. A user-provided scoring scheme overrides this.\n"
    "  options are:\n"
    "    -o <file_name>     output file [standard output, unless -bed is given]\n"
    "    -bed <file_name>   write BED6 to this file: target, start, end, query, per-mil identity, strand\n"
    "    -score <m>,<ti>,<tv>,<go>,<ge>  match, transition, transversion, gap open and gap extend scores,\n"
    "                       used everywhere, e.g. -score 1,-1,-2,-7,-1 [by default weak 2,-1,-3,-15,-1,\n"
    "                       medium 2,-2,-4,-20,-2 or strong 2,-3,-6,-25,-2, chosen by region composition]\n"
    "    -smin <n>          minimum score to report [50]\n"
    "    -xmin <n>          minimum x-drop: a match drops cells more than the x-drop below its best [200]\n"
    "    -xfac <x>          x-drop is the larger of xmin and xfac times the unit length [2.0]\n"
    "    -T <n>             number of scanning threads [8]\n" ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit(1) ; }

  FILE *outFile = 0, *bedFile = 0 ;	// outFile defaults to stdout unless -bed is given
  ScoreParams userScheme ;
  bool isUserScheme = false ;
  U64 scoreMin = 50 ;
  int xMin = 200 ;
  double xFac = 2.0 ;
  int nThread = 8 ;

  /* Process command line options */
  while (argc > 1 && **argv == '-')
    { if (!strcmp (*argv, "-o"))
        { if (!(outFile = fopen (*++argv, "w")))
            { fprintf (stderr, "failed to open output file %s\n", *argv) ;
              exit (1) ;
            }
          argc-- ;
        }
      else if (!strcmp (*argv, "-bed") && argc > 1)
        { if (!(bedFile = fopen (*++argv, "w")))
            { fprintf (stderr, "failed to open bed file %s\n", *argv) ;
              exit (1) ;
            }
          argc-- ;
        }
      else if (!strcmp (*argv, "-score") && argc > 1)
        { int m, ti, tv, go, ge ; char extra ;
          if (sscanf (*++argv, "%d,%d,%d,%d,%d%c", &m, &ti, &tv, &go, &ge, &extra) != 5)
            die ("bad -score %s - need 5 comma-separated integers, e.g. 1,-1,-2,-7,-1\n%s", *argv, usage) ;
          if (m <= 0) die ("bad match score %d in -score - must be +ve\n%s", m, usage) ;
          if (ti >= 0 || tv >= 0 || go >= 0 || ge >= 0)
            die ("bad -score %s - mismatch and gap scores must be -ve\n%s", *argv, usage) ;
          userScheme.match = m ; userScheme.transition = ti ; userScheme.transversion = tv ;
          userScheme.gapOpen = (I64)go * ONE ; userScheme.gapExtend = (I64)ge * ONE ;
          isUserScheme = true ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-T") && argc > 1)
        { nThread = atoi (*++argv) ;
          if (nThread < 1) die ("bad number of threads %s - must be at least 1\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-xmin") && argc > 1)
        { xMin = atoi (*++argv) ;
          if (xMin <= 0) die ("bad xmin %s - must be +ve\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-xfac") && argc > 1)
        { xFac = atof (*++argv) ;
          if (xFac < 0) die ("bad xfac %s - must be non-negative\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-smin") && argc > 1)
        { long long v = atoll (*++argv) ;
          if (v < 0) die ("bad smin score %s - must be non-negative\n%s", *argv, usage) ;
          scoreMin = v ;
          argc-- ;
        }
      else die ("unknown option: %s\n%s", *argv, usage) ;
      argc-- ; argv++ ;
    }
  if (argc != 2) die (usage) ;
  if (!outFile && !bedFile) outFile = stdout ;

  ScanParams scan ;
  scan.nScheme = isUserScheme ? 1 : MAX_SCHEME ;
  for (int k = 0 ; k < scan.nScheme ; ++k)
    { ScoreParams *sp = &scan.sp[k] ;
      *sp = isUserScheme ? userScheme : defaultScheme[k] ;
      for (int i = 4 ; i-- ;) for (int j = 4 ; j-- ;) baseScore[k][i][j] = (I64)sp->transversion * ONE ;
      baseScore[k][0][0] = baseScore[k][1][1] = baseScore[k][2][2] = baseScore[k][3][3] = (I64)sp->match * ONE + 1 ;
      baseScore[k][0][2] = baseScore[k][2][0] = baseScore[k][1][3] = baseScore[k][3][1] = (I64)sp->transition * ONE ;
    }
  if (isUserScheme) schemeChar[0] = 'u' ;
 
  Array aQuery = arrayCreate (256, Query) ;
  int nEntry = 0, nSkip = 0 ;
  OneSchema *schema = oneSchemaCreateFromText(repeatSchemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "lib", 1) ;
  if (of)
    { oneGoto (of, 'S', 1) ; oneReadLine (of) ; // start at first sequence
      while (of->lineType == 'S')
	{ Query *q = arrayp (aQuery, arrayMax(aQuery), Query) ;
	  bool isTandem = true ;
	  q->index = ++nEntry ;
	  q->len = oneLen(of) ;
	  q->seq = new (q->len, char) ;
	  char *s = oneDNAchar(of) ;
	  for (int i = 0 ; i < q->len ; ++i) q->seq[i] = dna2index4Conv[(int)s[i]] ;
	  while (oneReadLine (of) && of->lineType != 'S')
	    if (of->lineType == 'I' && !q->name)
	      { q->name = new (oneLen(of)+1, char) ; strcpy (q->name, oneString(of)) ; }
	    else if (of->lineType == 'Y' && strncmp (oneString(of), "tandem", 6))
	      isTandem = false ;
	  if (isTandem)
	    { if (!q->name)
		{ char buf[32] ; sprintf (buf, "tandem-%d", nEntry) ;
		  q->name = new (strlen(buf)+1, char) ; strcpy (q->name, buf) ;
		}
	    }
	  else
	    { ++nSkip ; newFree (q->seq, q->len, char) ;
	      if (q->name) { newFree (q->name, strlen(q->name)+1, char) ; q->name = 0 ; }
	      --arrayMax(aQuery) ; // fail to register this entry
	    }
	}
      oneFileClose (of) ;
    }
  else
    { SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
      if (!sio) die ("failed to open query sequence file %s", *argv) ;
      while (seqIOread (sio))
	{ Query *q = arrayp (aQuery, arrayMax(aQuery), Query) ;
	  q->len = sio->seqLen ; q->index = ++nEntry ;
	  q->name = new (strlen(sqioId(sio))+1, char) ; strcpy (q->name, sqioId(sio)) ;
	  q->seq = new (q->len, char) ; memcpy (q->seq, sqioSeq(sio), q->len) ;
	}
      seqIOclose (sio) ;
    }
  int nQuery = nQueryFwd = arrayMax(aQuery) ;
  if (!nQuery) die ("no tandem query sequences in %s", *argv) ;
  U64 totSat = 0 ;
  for (int i = 0 ; i < nQuery ; ++i) // reverse complements go after the forward queries
    { Query *qr = arrayp (aQuery, nQuery+i, Query) ;
      Query *q = arrp (aQuery, i, Query) ; // after arrayp(), which may have reallocated
      *qr = *q ; qr->isRC = true ;
      qr->seq = seqRevComp (q->seq, q->len) ;
      makeQueryScores (q, scan.nScheme) ; makeQueryScores (qr, scan.nScheme) ;
      I64 x = (I64)(xFac * q->len + 0.5) ;
      q->xDrop = qr->xDrop = (x > xMin ? x : xMin) * ONE ;
      totSat += q->len ;
    }
  --argc ; ++argv ;
  fprintf (stderr, "read in %d query sequences, total length %llu bp", nQuery, totSat) ;
  if (nSkip) fprintf (stderr, ", skipping %d whose type is not tandem", nSkip) ;
  fprintf (stderr, "\n") ;

  SeqIO *sio ;
  sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
  if (!sio) die ("failed to open target sequence file %s", *argv) ;
  scan.aQuery = aQuery ; scan.scoreMin = scoreMin ;
  Worker **worker = new0 (nThread, Worker*) ;
  for (int t = 0 ; t < nThread ; ++t) worker[t] = workerCreate (&scan) ;
  Array aTargetName = arrayCreate (256, char*) ;

  // with one thread the master scans directly, otherwise it feeds a queue serviced by nThread workers
  JobQueue queue ;
  ThreadArg *ta = 0 ;
  if (nThread > 1)
    { queueInit (&queue, nThread) ;
      ta = new0 (nThread, ThreadArg) ;
      for (int t = 0 ; t < nThread ; ++t)
	{ ta[t].scan = &scan ; ta[t].queue = &queue ; ta[t].worker = worker[t] ;
	  if (pthread_create (&ta[t].thread, 0, workerMain, &ta[t])) die ("failed to create thread %d", t) ;
	}
    }

  U64 nSeq = 0, totSeq = 0 ;
  while (seqIOread (sio))
    { ++nSeq ; totSeq += sio->seqLen ;
      if (sio->seqLen < 32) continue ;
      Job job ;
      job.target = arrayMax(aTargetName) ;
      job.name = new (strlen(sqioId(sio))+1, char) ; strcpy (job.name, sqioId(sio)) ;
      array(aTargetName, job.target, char*) = job.name ;
      job.len = sio->seqLen ;
      if (nThread == 1)
	{ job.seq = sqioSeq(sio) ; job.ownsSeq = false ;
	  workerRun (&scan, worker[0], &job) ;
	}
      else
	{ job.seq = new (job.len,char) ;   // copy, because the next seqIOread overwrites the buffer
	  memcpy (job.seq, sqioSeq(sio), job.len) ; job.ownsSeq = true ;
	  queuePush (&queue, &job) ;
	}
    }
  seqIOclose (sio) ;

  if (nThread > 1)
    { queueClose (&queue) ;
      for (int t = 0 ; t < nThread ; ++t) pthread_join (ta[t].thread, 0) ;
      queueDestroy (&queue) ;
      newFree (ta, nThread, ThreadArg) ;
    }

  // gather the matches from all the workers into the first worker's array, then sort and write
  U64 totStart = 0, totMatch = 0, totMatchLen = 0, totSameEnd = 0 ;
  Array aMatch = worker[0]->aMatch ;
  for (int t = 0 ; t < nThread ; ++t)
    { Worker *w = worker[t] ;
      totStart += w->totStart ; totMatch += w->totMatch ; totMatchLen += w->totMatchLen ; totSameEnd += w->totSameEnd ;
      if (t)
	for (U64 i = 0 ; i < arrayMax(w->aMatch) ; ++i)
	  array(aMatch, arrayMax(aMatch), Match*) = arr(w->aMatch, i, Match*) ;
    }
  arraySort (aMatch, matchOrder) ;
  for (U64 i = 0 ; i < arrayMax(aMatch) ; ++i)
    { Match *m = arr(aMatch, i, Match*) ;
      char *tName = arr(aTargetName, m->target, char*) ;
      Query *q    = arrp(aQuery, m->query, Query) ;
      char  *qName = q->name ;
      U64    count = m->bestScore & 0xffffffff ;
      char   strand = q->isRC ? '-' : '+' ;
      if (outFile)
	{ fprintf (outFile, "%s\t%llu\t%llu\t%s\t%llu\t%llu\t%c\t", tName, m->start, m->bestEnd, qName,
		   (U64)(m->bestScore >> 32), count, strand) ;
	  if (!m->nParse) fputc ('.', outFile) ; // just one copy
	  for (int k = 0 ; k < m->nParse ; ++k) fprintf (outFile, k ? ",%u" : "%u", m->parse[k]) ;
	  fputc ('\t', outFile) ; // the scheme(s) used
	  for (int k = 0 ; k < scan.nScheme ; ++k) if (m->schemes & (1 << k)) fputc (schemeChar[k], outFile) ;
	  fputc ('\n', outFile) ;
	}
      if (bedFile)
	fprintf (bedFile, "%s\t%llu\t%llu\t%s\t%d\t%c\n", tName, m->start, m->bestEnd, qName,
		 (int)(1000 * count / (m->bestEnd - m->start)), strand) ;
      newFree (m->parse, m->nParse ? m->nParse : 1, U32) ;
      newFree (m, 1, Match) ;
    }

  fprintf (stderr, "processed %llu target sequences, total length %llu bp\n", nSeq, totSeq) ;
  fprintf (stderr, "found %llu starts, %llu matches, total match length %llu bp\n",
	   totStart, totMatch, totMatchLen) ;
  if (totSameEnd) fprintf (stderr, "(after removing %llu matches ending where a better one of the same query ends)\n", totSameEnd) ;
  if (scan.nScheme > 1) // regions scanned with a stricter scheme, or walled, for each query strand
    for (int j = 0 ; j < arrayMax(aQuery) ; ++j)
      { U64 n[MAX_SCHEME+1], len[MAX_SCHEME+1] ;
	for (int k = 0 ; k <= MAX_SCHEME ; ++k)
	  { n[k] = len[k] = 0 ;
	    for (int t = 0 ; t < nThread ; ++t)
	      { n[k] += worker[t]->regionN[j*(MAX_SCHEME+1)+k] ; len[k] += worker[t]->regionLen[j*(MAX_SCHEME+1)+k] ; }
	  }
	if (!n[1] && !n[2] && !n[3]) continue ;
	Query *q = arrp(aQuery, j, Query) ;
	fprintf (stderr, "%s %c: medium %llu regions %llu bp, strong %llu regions %llu bp, walled %llu regions %llu bp\n",
		 q->name, q->isRC ? '-' : '+', n[1], len[1], n[2], len[2], n[3], len[3]) ;
      }
  timeTotal(stderr) ;

  if (outFile && outFile != stdout) fclose (outFile) ;
  if (bedFile) fclose (bedFile) ;
  for (U64 i = 0 ; i < arrayMax(aTargetName) ; ++i)
    { char *name = arr(aTargetName, i, char*) ; newFree (name, strlen(name)+1, char) ; }
  arrayDestroy (aTargetName) ;
  for (int t = 0 ; t < nThread ; ++t) workerDestroy (worker[t], &scan) ;
  newFree (worker, nThread, Worker*) ;
  queryDestroy (aQuery) ;
  oneSchemaDestroy (schema) ;
}
