/*  File: satmatch.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description:
 * Pairwise DNA sequence alignment - a simple implementation
 *
 * Usage: satmatch <seqfile> [-o outfile] [-u unit] [-s seqfile] [-c coord]
 *
 * Exports functions:
 *   - loadSequences(): loads all sequences from a file
 *   - pairwiseAlign(): performs pairwise alignment between two sequences
 *   - satmatch(): main alignment routine
 *
 * HISTORY:
 * Last edited: Oct 27 01:34 2026 (claude)
 * Created: Tue Aug 11 08:32:53 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h" // for VERSION
#include "onerep.h"   // for the .1lib schema
#include "satscore.h" // satscan's scoring, used to link representatives
#include "array.h"
#include "seqio.h"
#include <string.h>
#include <pthread.h>

typedef struct {
  int   len ;
  char *name ;
  char *seq ;
  char *rSeq ;
  int   level ;	// the least strict satscan scheme safe for its composition, MAX_SCHEME if none is
} Seq ;

typedef struct {
  int transition ;      /* A<->G, C<->T */
  int transversion ;    /* A<->C, A<->T, G<->C, G<->T */
  int gap ;             /* gap penalty per position */
} ScoreParams ;

static inline int min (int x, int y) { return x < y ? x : y ; }

#define NO_THRESHOLD (1 << 30)

// The best (lowest) cost of aligning all of the longer sequence against the shorter one as a circle. All
// costs are non-negative, so the lowest cost in a row can only rise in later rows: once it is above
// threshold the alignment stops, returning that lowest cost, a lower bound above threshold, not the score.
int needlemanWunsch (Seq *s1, Seq *s2, ScoreParams *sp, int threshold, int *jMax)
{
  if (s2->len > s1->len) { Seq *s = s1 ; s1 = s2 ; s2 = s ; } // swap
  int len1 = s1->len ;
  int len2 = s2->len ;
  
  // Build 4x4 cost table: index encoding a=0,c=1,g=2,t=3
  int cost[4][4] ;
  int i, j ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      cost[i][j] = sp->transversion ;
  cost[0][0] = cost[1][1] = cost[2][2] = cost[3][3] = 0 ;  // Diagonal = match
  cost[0][2] = cost[2][0] = cost[1][3] = cost[3][1] = sp->transition ;  // A<>G (0^2), C<>T (1^3)

  // Allocate two rows; the first is 0, so the alignment can start at any rotation of s2
  int *prev = new0 (len2+1, int) ;
  int *row = new (len2+1, int) ;

  // Fill DP
  //  for (j = 0 ; j < len1 ; ++j) printf ("  %c", s2->seq[j]) ; putchar ('\n') ;
  for (i = 1 ; i <= len1 ; i++)
    { // s2 is circular: column 0 continues from the end of s2, so the shorter s2 can wrap many times
      row[0] = min (prev[len2-1]+cost[s1->seq[i-1]][s2->seq[len2-1]], sp->gap+prev[len2]) ;
      int rowMin = row[0] ;
      // printf ("%c %2d", index2char[s1->seq[i-1]], row[0]) ;
      for (j = 1 ; j <= len2 ; j++)
        { int d = prev[j-1] + cost[s1->seq[i-1]][s2->seq[j-1]] ;  // diagonal: match/mismatch
          int g = sp->gap + min (prev[j], row[j-1]) ;             // up or left + gap
          row[j] = min (d, g) ;
	  if (row[j] < rowMin) rowMin = row[j] ;
	  // printf (" %2d", row[j]) ;
        }
      // putchar ('\n') ;
      int *tmp = prev ; prev = row ; row = tmp ; // swap rows
      if (rowMin > threshold) // no alignment can now end at or below threshold
	{ newFree (prev, len2+1, int) ; newFree (row, len2+1, int) ;
	  return rowMin ;
	}
    }
  
  int score = 3 * len1 ; // the best end, at any position in s2
  for (j = 0 ; j < len2 ; ++j)
    if (prev[j] < score) { score = prev[j] ; if (jMax) *jMax = j ; }
  
  newFree (prev, len2+1, int) ;
  newFree (row, len2+1, int) ;

  return score ;
}

// Single-cycle alignment, used by the greedy clustering: the best cost of aligning all of the shorter
// sequence S against exactly one cycle of the longer one R, starting at any rotation of R, so with no
// multiple copies. R is on the columns, as a circle, and S on the rows. Each cell holds the cost in its upper
// 32 bits and, in its lower 32 bits, the column of R at which its least-cost path started, plus L each time
// the path crossed from the last column to the first, carried along by taking the minimum of the packed
// values. A path starting at column 0 begins from the cell for column L-1, so every single cycle crosses
// exactly once. On the last row, a cell whose start is the column after it, having crossed once, has gone
// exactly once round R (one that has not crossed has consumed no column of R), so its cost is that of a real single-cycle alignment: the least of these is an
// upper bound UB on the best. The least cost of the whole row, ending anywhere, is a lower bound LB, as the
// paths it allows include all single-cycle ones. So if LB > threshold the pair cannot link, and if
// UB <= threshold it can. Otherwise, rarely, R is aligned globally from each of the starts of the cheapest
// last-row cells. All costs are non-negative, so once a row's least cost is above threshold the alignment
// stops, returning it, a bound above threshold.

#define NTRY 64 // the most last-row cells whose starts are tried when the bounds leave a pair undecided
typedef struct { I64 aligned, earlyExit, lowerBound, upperBound, ambiguous, resolved ; } CyclicStats ; // counts of alignments, kept atomically
static CyclicStats cyclicStats ;

static void costTable (ScoreParams *sp, int cost[4][4])
{
  for (int i = 0 ; i < 4 ; i++) for (int j = 0 ; j < 4 ; j++) cost[i][j] = sp->transversion ;
  cost[0][0] = cost[1][1] = cost[2][2] = cost[3][3] = 0 ;
  cost[0][2] = cost[2][0] = cost[1][3] = cost[3][1] = sp->transition ;
}

// an ordinary global alignment of all of S against all of R, starting at column start of R
static int rotatedCost (char *s, int n, char *r, int L, int start, ScoreParams *sp, int threshold)
{
  int cost[4][4] ; costTable (sp, cost) ;
  int *prev = new (L+1, int), *cur = new (L+1, int) ;
  for (int j = 0 ; j <= L ; ++j) prev[j] = j * sp->gap ;
  int result = 0 ;
  for (int i = 0 ; i < n ; ++i)
    { cur[0] = (i+1) * sp->gap ;
      int rowMin = cur[0] ;
      for (int j = 1, k = start ; j <= L ; ++j, k = (k+1 == L) ? 0 : k+1)
	{ int d = prev[j-1] + cost[(int)s[i]][(int)r[k]], g = sp->gap + min (prev[j], cur[j-1]) ;
	  cur[j] = min (d, g) ;
	  if (cur[j] < rowMin) rowMin = cur[j] ;
	}
      int *t = prev ; prev = cur ; cur = t ;
      if (rowMin > threshold) { result = rowMin ; break ; }
      result = prev[L] ;
    }
  newFree (prev, L+1, int) ; newFree (cur, L+1, int) ;
  return result ;
}

static int cyclicCost (Seq *a, Seq *b, ScoreParams *sp, int threshold)
{
  Seq *rs = a->len >= b->len ? a : b, *ss = a->len >= b->len ? b : a ;
  char *r = rs->seq, *s = ss->seq ;
  int L = rs->len, n = ss->len ;
  int cost[4][4] ; costTable (sp, cost) ;
  I64 G = (I64)sp->gap << 32 ;
  I64 WRAP = L ; // added to the start label on crossing from column L-1 to 0, so labels count cycles
  I64 *prev = new (L, I64), *cur = new (L, I64) ;
  for (int j = 0 ; j < L ; ++j) prev[j] = (j+1 == L) ? 0 : j+1 ; // cost 0, starting at the next column
  __atomic_fetch_add (&cyclicStats.aligned, 1, __ATOMIC_RELAXED) ;
  for (int i = 0 ; i < n ; ++i)
    { int c = s[i] ;
      I64 rowMin ;
      for (int j = 0 ; j < L ; ++j) // diagonal from column j-1, and a gap in R (S[i] alone) from above
	{ I64 d = (j ? prev[j-1] : prev[L-1] + WRAP) + ((I64)cost[c][(int)r[j]] << 32), v = prev[j] + G ;
	  cur[j] = d < v ? d : v ;
	}
      // gaps in S along the row (R[j] alone), round the circle: a first pass from column 0, entering from
      // column L-1, then on round while it still improves cells
      I64 h = cur[L-1] + G + WRAP ;
      for (int j = 0 ; j < L ; ++j) { if (h < cur[j]) cur[j] = h ; h = cur[j] + G ; }
      h += WRAP ; // round again
      for (int j = 0 ; j < L && h < cur[j] ; ++j) { cur[j] = h ; h += G ; }
      rowMin = cur[0] ;
      for (int j = 1 ; j < L ; ++j) if (cur[j] < rowMin) rowMin = cur[j] ;
      I64 *t = prev ; prev = cur ; cur = t ;
      if ((rowMin >> 32) > threshold) // LB > threshold: no alignment can now end at or below it
	{ __atomic_fetch_add (&cyclicStats.earlyExit, 1, __ATOMIC_RELAXED) ;
	  newFree (prev, L, I64) ; newFree (cur, L, I64) ;
	  return rowMin >> 32 ;
	}
    }
  int lb = NO_THRESHOLD, ub = NO_THRESHOLD ;
  for (int j = 0 ; j < L ; ++j)
    { int c = prev[j] >> 32, label = prev[j] & 0xffffffff, start = label % L, wraps = label / L ;
      if (c < lb) lb = c ;
      if (start == ((j+1 == L) ? 0 : j+1) && wraps == 1 && c < ub) ub = c ; // exactly once round R
    }
  int result = ub ;
  if (lb > threshold) __atomic_fetch_add (&cyclicStats.lowerBound, 1, __ATOMIC_RELAXED) ;
  else if (ub <= threshold) __atomic_fetch_add (&cyclicStats.upperBound, 1, __ATOMIC_RELAXED) ;
  if (lb <= threshold && ub > threshold) // undecided: align from the starts suggested by the cheapest cells
    { __atomic_fetch_add (&cyclicStats.ambiguous, 1, __ATOMIC_RELAXED) ;
      // Paths from neighbouring starts merge within a few rows, so a cell often holds a neighbour's start.
      // For each of the cheapest cells, at column j, in order, try the start j+1 that would complete a cycle
      // there, and the start it holds, until the pair links or the cells are above threshold: the cost of a
      // cycle ending at a cell is at least the cell's.
      int tried[2*NTRY], nTried = 0 ;
      for (int k = 0 ; k < NTRY && result > threshold ; ++k)
	{ int bj = -1 ; // the cheapest cell not yet used
	  for (int j = 0 ; j < L ; ++j) if ((prev[j] >> 32) < NO_THRESHOLD && (bj < 0 || prev[j] < prev[bj])) bj = j ;
	  if (bj < 0 || (prev[bj] >> 32) > threshold) break ;
	  int starts[2] = { (bj+1) % L, (prev[bj] & 0xffffffff) % L } ;
	  prev[bj] = (I64)NO_THRESHOLD << 32 ; // used
	  for (int q = 0 ; q < 2 && result > threshold ; ++q)
	    { int t ; for (t = 0 ; t < nTried ; ++t) if (tried[t] == starts[q]) break ;
	      if (t < nTried) continue ;
	      tried[nTried++] = starts[q] ;
	      int x = rotatedCost (s, n, r, L, starts[q], sp, threshold) ;
	      if (x < result) result = x ;
	    }
	}
      if (result <= threshold) __atomic_fetch_add (&cyclicStats.resolved, 1, __ATOMIC_RELAXED) ;
    }
  newFree (prev, L, I64) ; newFree (cur, L, I64) ;
  return result ;
}

Array randomSeqs (int n, int len)
{
  Array a = arrayCreate (n, Seq) ;
  for (int i = 0 ; i < n ; ++i)
    { Seq *s = arrayp (a, arrayMax(a), Seq) ;
      s->len = len ;
      s->name = new(32,char) ; sprintf (s->name, "random%d-%d", len, i+1) ;
      s->seq = new (s->len, char) ; for (int j = 0 ; j < len ; ++j) s->seq[j] = rand()&0x3 ;
      s->rSeq = seqRevComp (s->seq, s->len) ;
      s->level = satSchemeLevel (s->seq, s->len) ;
    }
  return a ;
}

// All pairs are scored by nThread threads, each taking the next unscored row i and scoring (i,j) for all
// j > i, before the pairs are linked into clusters in order, so the results do not depend on nThread.
typedef struct {
  Array aSeq ;
  ScoreParams *sp ;
  int nSeq ;
  int *scoreij ;	// nSeq*nSeq: the better score of (i,j) over both orientations of j, for j > i
  int *offsetij ;	// if not 0, the end offset of that alignment, as -(offset+1) if j was reverse complemented
  double threshFac ;	// if > 0, alignments stop once they cannot link, i.e. are above pairThreshold()
  int nextRow ;
  pthread_mutex_t lock ;
} PairJob ;

// two sequences are linked if their score is at most this
static inline int pairThreshold (double threshFac, Seq *a, Seq *b)
{ return threshFac * (a->len < b->len ? b->len : a->len) - 10 ; }

// the better score of i against j in either orientation, or a bound above threshold if both are above it;
// *offset as for offsetij. The second orientation only needs to beat the first.
static int pairScore (Array aSeq, ScoreParams *sp, int i, int j, int threshold, int *offset)
{
  Seq *si = arrp(aSeq,i,Seq), *sj = arrp(aSeq,j,Seq) ;
  Seq sjRC = *sj ; sjRC.seq = sj->rSeq ; // so sj itself is never changed while other threads read it
  int jMax1 = 0, jMax2 = 0 ;
  int score1 = needlemanWunsch (si, sj, sp, threshold, &jMax1) ;
  int score2 = needlemanWunsch (si, &sjRC, sp, min (threshold, score1), &jMax2) ;
  if (offset) *offset = (score1 <= score2) ? jMax1 : -(jMax2+1) ;
  return min (score1, score2) ;
}

// the better single-cycle cost of i against j in either orientation, or a bound above threshold
static int cyclicPairScore (Array aSeq, ScoreParams *sp, int i, int j, int threshold)
{
  Seq *si = arrp(aSeq,i,Seq), *sj = arrp(aSeq,j,Seq) ;
  Seq sjRC = *sj ; sjRC.seq = sj->rSeq ;
  int score1 = cyclicCost (si, sj, sp, threshold) ;
  int score2 = cyclicCost (si, &sjRC, sp, min (threshold, score1)) ;
  return min (score1, score2) ;
}

static void scorePair (PairJob *pj, int i, int j)
{
  size_t k = (size_t)i*pj->nSeq + j ;
  int threshold = pj->threshFac > 0 ?
    pairThreshold (pj->threshFac, arrp(pj->aSeq,i,Seq), arrp(pj->aSeq,j,Seq)) : NO_THRESHOLD ;
  pj->scoreij[k] = pairScore (pj->aSeq, pj->sp, i, j, threshold, pj->offsetij ? &pj->offsetij[k] : 0) ;
}

static void *pairThread (void *arg)
{
  PairJob *pj = (PairJob*) arg ;
  while (true)
    { pthread_mutex_lock (&pj->lock) ;
      int i = pj->nextRow++ ;
      pthread_mutex_unlock (&pj->lock) ;
      if (i >= pj->nSeq) break ;
      for (int j = i+1 ; j < pj->nSeq ; ++j) scorePair (pj, i, j) ;
    }
  return 0 ;
}

static void runThreads (void *(*f)(void*), void *arg, int nThread)
{
  if (nThread == 1) { f (arg) ; return ; }
  pthread_t *thread = new (nThread, pthread_t) ;
  for (int t = 0 ; t < nThread ; ++t)
    if (pthread_create (&thread[t], 0, f, arg)) die ("failed to create thread %d", t) ;
  for (int t = 0 ; t < nThread ; ++t) pthread_join (thread[t], 0) ;
  newFree (thread, nThread, pthread_t) ;
}

// a list of pairs to score in parallel, each with i < j as in the all-pairs scoring; only pairs whose
// score is UNSCORED are scored, each with its own threshold: with satscan's scheme level if level >= 0,
// by satPairScore(), else with satmatch's costs, by single-cycle alignment if level is CYCLIC, else
// with the shorter sequence wrapping round any number of times
#define UNSCORED (-(1 << 30))
#define CYCLIC (-2)
typedef struct { int i, j, threshold, score, offset, level ; } Pair ;

// the best satscan score of a single-cycle alignment of i and j, from any rotation of the longer but with no
// multiple copies, in either orientation, using scheme level; a bound at most threshold if it cannot be above
static int satPairScore (Array aSeq, int i, int j, int level, int threshold)
{
  Seq *a = arrp(aSeq,i,Seq), *b = arrp(aSeq,j,Seq) ;
  const SatScheme *sp = &satScheme[level] ;
  int score1 = satCyclicScore (a->seq, a->len, b->seq, b->len, sp, threshold) ;
  int score2 = satCyclicScore (a->seq, a->len, b->rSeq, b->len, sp, score1 > threshold ? score1 : threshold) ;
  return score1 > score2 ? score1 : score2 ;
}

typedef struct {
  Array aSeq ;
  ScoreParams *sp ;
  Pair *pair ;
  I64 nPair, next ;
  pthread_mutex_t lock ;
} PairListJob ;

static void *pairListThread (void *arg)
{
  PairListJob *pl = (PairListJob*) arg ;
  static const int CHUNK = 16 ;
  while (true)
    { pthread_mutex_lock (&pl->lock) ;
      I64 k0 = pl->next ; pl->next += CHUNK ;
      pthread_mutex_unlock (&pl->lock) ;
      if (k0 >= pl->nPair) break ;
      for (I64 k = k0 ; k < k0 + CHUNK && k < pl->nPair ; ++k)
	{ Pair *p = &pl->pair[k] ;
	  if (p->score == UNSCORED)
	    p->score = p->level >= 0 ? satPairScore (pl->aSeq, p->i, p->j, p->level, p->threshold)
	      : p->level == CYCLIC ? cyclicPairScore (pl->aSeq, pl->sp, p->i, p->j, p->threshold)
	      : pairScore (pl->aSeq, pl->sp, p->i, p->j, p->threshold, &p->offset) ;
	}
    }
  return 0 ;
}

static void scorePairs (Array aSeq, ScoreParams *sp, Array aPair, int nThread)
{
  PairListJob pl = { aSeq, sp, arrp(aPair,0,Pair), arrayMax(aPair), 0 } ;
  pthread_mutex_init (&pl.lock, 0) ;
  runThreads (pairListThread, &pl, nThread) ;
  pthread_mutex_destroy (&pl.lock) ;
}

static inline Pair *addPair (Array aPair, int i, int j, int threshold)
{ Pair *p = arrayp(aPair, arrayMax(aPair), Pair) ;
  if (i < j) { p->i = i ; p->j = j ; } else { p->i = j ; p->j = i ; }
  p->threshold = threshold ; p->score = UNSCORED ; p->offset = 0 ; p->level = -1 ;
  return p ;
}

// A pair for the greedy clustering, scored by single-cycle alignment unless its lengths alone show that it
// cannot link: every alignment has at least |lengthDifference| gap positions, and all costs are
// non-negative, so its cost is at least gap * |lengthDifference|.
static I64 nLengthReject = 0 ;
static void greedyPair (Array aSeq, ScoreParams *sp, Array aPair, int i, int j, int threshold)
{
  Pair *p = addPair (aPair, i, j, threshold) ;
  int bound = sp->gap * abs (arrp(aSeq,i,Seq)->len - arrp(aSeq,j,Seq)->len) ;
  if (bound > threshold) { p->score = bound ; ++nLengthReject ; } // not aligned
  else p->level = CYCLIC ;
}

static Array lengthOrderSeqs ; // for lengthOrder(), which sorts sequence indices by length, then index
static int lengthOrder (const void *a, const void *b)
{ int x = *(const int*)a, y = *(const int*)b ;
  int lx = arrp(lengthOrderSeqs,x,Seq)->len, ly = arrp(lengthOrderSeqs,y,Seq)->len ;
  return lx != ly ? (lx < ly ? -1 : 1) : (x < y ? -1 : x > y) ;
}

// Greedy clustering: the sequences are taken from shortest to longest, and each joins the
// cluster of the representative it scores best against by single-cycle alignment, cyclicCost(), if that
// is within the threshold for its length,
// else becomes the representative of a new cluster. So each representative is the shortest member of its
// cluster. To use threads, a batch of sequences is scored against all the representatives found before
// the batch and against all the earlier sequences in the batch, and the batch is then resolved in order,
// which gives the same clusters as taking the sequences one at a time. These representatives have all
// been compared with each other already, so linkRepresentatives() uses the ones chosen at the end.
static int clusterGreedy (Array aSeq, ScoreParams *sp, double threshFac, int nThread, int *cluster)
{
  int nSeq = arrayMax(aSeq) ;
  int *order = new (nSeq, int) ;
  for (int i = 0 ; i < nSeq ; ++i) order[i] = i ;
  lengthOrderSeqs = aSeq ; qsort (order, nSeq, sizeof(int), lengthOrder) ;
  Array aRep = arrayCreate (1024, int) ; // representatives, in order of cluster number - 1
  Array aPair = arrayCreate (1024, Pair) ;
  int batch = 32 * nThread ;
  for (int b0 = 0, b1 ; b0 < nSeq ; b0 = b1)
    { b1 = b0 + batch ; if (b1 > nSeq) b1 = nSeq ;
      int nRep0 = arrayMax(aRep) ;
      arrayMax(aPair) = 0 ;
      for (int k = b0 ; k < b1 ; ++k)
	{ Seq *s = arrp(aSeq,order[k],Seq) ;
	  int threshold = threshFac * s->len - 10 ; // s is at least as long as any representative
	  for (int r = 0 ; r < nRep0 ; ++r) greedyPair (aSeq, sp, aPair, order[k], arr(aRep,r,int), threshold) ;
	  for (int k2 = b0 ; k2 < k ; ++k2) greedyPair (aSeq, sp, aPair, order[k], order[k2], threshold) ;
	}
      scorePairs (aSeq, sp, aPair, nThread) ;
      Pair *p = arrp(aPair,0,Pair) ;
      for (int k = b0 ; k < b1 ; ++k)
	{ int i = order[k] ; Seq *s = arrp(aSeq,i,Seq) ;
	  int threshold = p->threshold ; // the same for all of i's pairs
	  int best = 1 << 30, bestCluster = 0, bestRep = -1 ;
	  for (int r = 0 ; r < nRep0 ; ++r, ++p)
	    if (p->score <= threshold && p->score < best)
	      { best = p->score ; bestCluster = r+1 ; bestRep = arr(aRep,r,int) ; }
	  for (int k2 = b0 ; k2 < k ; ++k2, ++p) // later clusters, so a tie keeps the earlier one
	    if (cluster[order[k2]] > nRep0 && arr(aRep,cluster[order[k2]]-1,int) == order[k2]
		&& p->score <= threshold && p->score < best)
	      { best = p->score ; bestCluster = cluster[order[k2]] ; bestRep = order[k2] ; }
	  if (bestCluster)
	    { cluster[i] = bestCluster ;
	      printf ("  link cluster %d: %d (%d) to representative %d (%d) score %d thresh %d\n",
		      bestCluster, i, s->len, bestRep, arrp(aSeq,bestRep,Seq)->len, best, threshold) ;
	    }
	  else
	    { array(aRep, arrayMax(aRep), int) = i ;
	      cluster[i] = arrayMax(aRep) ;
	      printf ("  new cluster %d: %d (%d)\n", cluster[i], i, s->len) ;
	    }
	}
    }
  int nCluster = arrayMax(aRep) ;
  fprintf (stderr, "greedy clustering: %d clusters; %lld pairs rejected by length, %lld alignments: rejected early %lld, by lower bound %lld, accepted by upper bound %lld, undecided %lld of which %lld linked by rotated alignment\t",
	   nCluster, (long long)nLengthReject, (long long)cyclicStats.aligned, (long long)cyclicStats.earlyExit, (long long)cyclicStats.lowerBound,
	   (long long)cyclicStats.upperBound, (long long)cyclicStats.ambiguous, (long long)cyclicStats.resolved) ;
  timeUpdate (stderr) ;

  arrayDestroy (aPair) ; arrayDestroy (aRep) ; newFree (order, nSeq, int) ;
  return nCluster ;
}

// The representative of each cluster is the member with the lowest total score to all the members,
// among those at most 1.1 times as long as the shortest member, the candidates. (Taking the lowest
// overall, and only if that is too long the lowest among the candidates, gives the same choice.) So
// only the candidates need scoring against all the members. This chooses rep[k] for each cluster k
// with isDo[k] (all if isDo is 0), and sets minLen[k] and, for its members, repScore[i], the score of
// i to rep[k]. With -allPairs, scoreij has the scores of all pairs, exact if isExact, else exact only
// within the threshold.
static void chooseRepresentatives (Array aSeq, ScoreParams *sp, int nThread, int *cluster, int nCluster,
				   bool *isDo, int *scoreij, bool isExact, double threshFac,
				   int *minLen, int *rep, int *repScore)
{
  int nSeq = arrayMax(aSeq) ;
  for (int k = 0 ; k <= nCluster ; ++k) if (!isDo || isDo[k]) { minLen[k] = 1 << 30 ; rep[k] = -1 ; }
  for (int i = 0 ; i < nSeq ; ++i)
    if (cluster[i] && (!isDo || isDo[cluster[i]]) && arrp(aSeq,i,Seq)->len < minLen[cluster[i]])
      minLen[cluster[i]] = arrp(aSeq,i,Seq)->len ;
  bool *isCand = new0 (nSeq, bool) ;
  Array *member = new0 (nCluster+1, Array) ; // the members of each cluster, in input order
  for (int i = 0 ; i < nSeq ; ++i)
    if (cluster[i] && (!isDo || isDo[cluster[i]]))
      { if (arrp(aSeq,i,Seq)->len <= 1.1*minLen[cluster[i]]) isCand[i] = true ;
	if (!member[cluster[i]]) member[cluster[i]] = arrayCreate (8, int) ;
	array(member[cluster[i]], arrayMax(member[cluster[i]]), int) = i ;
      }
  Array aPair = arrayCreate (1024, Pair) ;
  for (int k = 1 ; k <= nCluster ; ++k)
    if (member[k])
      for (int a = 0 ; a < arrayMax(member[k]) ; ++a)
	for (int b = a+1 ; b < arrayMax(member[k]) ; ++b)
	  { int i = arr(member[k],a,int), j = arr(member[k],b,int) ;
	    if (!isCand[i] && !isCand[j]) continue ;
	    Pair *p = addPair (aPair, i, j, NO_THRESHOLD) ;
	    if (scoreij)
	      { int score = scoreij[(size_t)p->i*nSeq+p->j] ;
		if (isExact || score <= pairThreshold (threshFac, arrp(aSeq,i,Seq), arrp(aSeq,j,Seq)))
		  p->score = score ; // else it may be a bound, so it is rescored
	      }
	  }
  scorePairs (aSeq, sp, aPair, nThread) ; // the pairs not yet scored exactly
  int *sumScore = new0 (nSeq, int) ;
  for (I64 k = 0 ; k < arrayMax(aPair) ; ++k)
    { Pair *p = arrp(aPair,k,Pair) ;
      if (isCand[p->i]) sumScore[p->i] += p->score ;
      if (isCand[p->j]) sumScore[p->j] += p->score ;
    }
  for (int k = 1 ; k <= nCluster ; ++k)
    if (member[k])
      for (int a = 0 ; a < arrayMax(member[k]) ; ++a)
	{ int i = arr(member[k],a,int) ;
	  if (isCand[i] && (rep[k] < 0 || sumScore[i] < sumScore[rep[k]])) rep[k] = i ;
	}
  for (I64 k = 0 ; k < arrayMax(aPair) ; ++k)
    { Pair *p = arrp(aPair,k,Pair) ;
      if (rep[cluster[p->i]] == p->i) repScore[p->j] = p->score ;
      else if (rep[cluster[p->j]] == p->j) repScore[p->i] = p->score ;
    }
  for (int k = 1 ; k <= nCluster ; ++k)
    if (member[k]) { repScore[rep[k]] = 0 ; arrayDestroy (member[k]) ; }
  fprintf (stderr, "choosing representatives: %lld pairs scored\t", (long long)arrayMax(aPair)) ; timeUpdate (stderr) ;
  newFree (member, nCluster+1, Array) ; newFree (isCand, nSeq, bool) ; newFree (sumScore, nSeq, int) ;
  arrayDestroy (aPair) ;
}

// Single linkage of the clusters over all pairs of their representatives, in cluster order, keeping the
// lower cluster number, so that e.g. subfamilies of a satellite are joined. Pairs are scored as satscan
// would, more sensitively than by satmatch's costs: a single-cycle alignment, from any rotation of the
// longer representative but with no multiple copies of the shorter, with the stricter of the two schemes
// safe for their compositions, and are linked if the score is above minScore. A pair whose lengths alone
// show it cannot score above minScore, satCyclicBound(), is not aligned. Representatives safe for no scheme are not linked. Sets isMerged[k] for
// each cluster k that others are merged into, and returns the number of merges.
static int linkRepresentatives (Array aSeq, ScoreParams *sp, int minScore, int nThread,
				int *cluster, int nCluster, int *rep, bool *isMerged)
{
  int nSeq = arrayMax(aSeq) ;
  Array aK = arrayCreate (nCluster, int) ; // the clusters with members
  bool *isUsed = new0 (nCluster+1, bool) ;
  for (int i = 0 ; i < nSeq ; ++i) if (cluster[i]) isUsed[cluster[i]] = true ;
  for (int k = 1 ; k <= nCluster ; ++k) if (isUsed[k]) array(aK, arrayMax(aK), int) = k ;
  int nK = arrayMax(aK) ;
  Array aPair = arrayCreate (1024, Pair) ;
  int nBound = 0 ;
  long long aligned0 = satCyclicAligned, undecided0 = satCyclicUndecided, resolved0 = satCyclicResolved ;
  long long early0 = satCyclicEarly, upper0 = satCyclicUpper, lower0 = satCyclicLower ;
  for (int a = 0 ; a < nK ; ++a)
    for (int b = a+1 ; b < nK ; ++b)
      { int i = rep[arr(aK,a,int)], j = rep[arr(aK,b,int)] ;
	int li = arrp(aSeq,i,Seq)->level, lj = arrp(aSeq,j,Seq)->level ;
	Pair *p = addPair (aPair, i, j, minScore) ;
	if (li == MAX_SCHEME || lj == MAX_SCHEME) p->score = minScore ; // not linked, and not aligned
	else
	  { int level = li > lj ? li : lj ;
	    int bound = satCyclicBound (arrp(aSeq,i,Seq)->len, arrp(aSeq,j,Seq)->len, &satScheme[level]) ;
	    if (bound <= minScore) { p->score = bound ; ++nBound ; } // can't link: not aligned
	    else p->level = level ;
	  }
      }
  scorePairs (aSeq, sp, aPair, nThread) ;
  int *map = new (nCluster+1, int) ; // the cluster each cluster is merged into
  for (int k = 0 ; k <= nCluster ; ++k) map[k] = k ;
  int nMerge = 0 ;
  Pair *p = arrp(aPair,0,Pair) ;
  for (int a = 0 ; a < nK ; ++a)
    for (int b = a+1 ; b < nK ; ++b, ++p)
      { int ka = arr(aK,a,int), kb = arr(aK,b,int) ;
	if (p->score > minScore && map[ka] != map[kb])
	  { int to = min (map[ka], map[kb]), from = map[ka] + map[kb] - to ;
	    printf ("  merge cluster %d to %d via representatives %d (%d) and %d (%d) score %d scheme %c minScore %d\n",
		    from, to, rep[ka], arrp(aSeq,rep[ka],Seq)->len, rep[kb], arrp(aSeq,rep[kb],Seq)->len,
		    p->score, satSchemeChar[p->level], minScore) ;
	    for (int k = 1 ; k <= nCluster ; ++k) if (map[k] == from) map[k] = to ;
	    isMerged[to] = true ;
	    ++nMerge ;
	  }
      }
  for (int i = 0 ; i < nSeq ; ++i) cluster[i] = map[cluster[i]] ;
  fprintf (stderr, "linking representatives: %d merges, %d clusters; %d pairs rejected by length, %lld alignments: rejected early %lld, by upper bound %lld, accepted by lower bound %lld, undecided %lld of which %lld linked by rotated alignment\t",
	   nMerge, nK - nMerge, nBound, satCyclicAligned - aligned0, satCyclicEarly - early0, satCyclicUpper - upper0, satCyclicLower - lower0,
	   satCyclicUndecided - undecided0, satCyclicResolved - resolved0) ;
  timeUpdate (stderr) ;
  newFree (map, nCluster+1, int) ; newFree (isUsed, nCluster+1, bool) ;
  arrayDestroy (aPair) ; arrayDestroy (aK) ;
  return nMerge ;
}

// The number of copies of the unit in the array a sequence came from, as int(array length / unit length),
// with the array read from a tancons name, cons<unit>#<seq>:<start>-<end>#<count>#<identity>; else 0.
static int arrayCopies (Seq *s)
{
  char *a = strchr (s->name, '#') ;
  if (!a || strncmp (s->name, "cons", 4)) return 0 ;
  int unit = atoi (s->name + 4) ;
  char *b = strchr (a+1, '#') ;
  if (!b || unit <= 0) return 0 ;
  char *c = b ; while (c > a && *c != ':') --c ; // the last ':' before the second '#', as seq may contain ':'
  long long start, end ;
  if (*c != ':' || sscanf (c+1, "%lld-%lld", &start, &end) != 2 || end <= start) return 0 ;
  return (end - start) / unit ;
}

static int intDecreasing (const void *a, const void *b)
{ int x = *(const int*)a, y = *(const int*)b ; return x > y ? -1 : x < y ; }

// the N50 of the array copy numbers of the members of a cluster: the copy number at which the members
// taken from the most copies down first include half of all their copies
static int n50Copies (Array aSeq, Array members)
{
  int n = arrayMax(members), *c = new (n, int) ;
  long long total = 0, sum = 0 ;
  for (int a = 0 ; a < n ; ++a) total += (c[a] = arrayCopies (arrp(aSeq, arr(members,a,int), Seq))) ;
  qsort (c, n, sizeof(int), intDecreasing) ;
  int x = 0 ;
  for (int a = 0 ; a < n && 2*sum < total ; ++a) { sum += c[a] ; x = c[a] ; }
  newFree (c, n, int) ;
  return x ;
}

void seqArrayDestroy (Array a)
{
  for (int i = 0 ; i < arrayMax(a) ; ++i)
    { Seq *s = arrp(a,i,Seq) ;
      newFree (s->seq, s->len, char) ;
      if (s->rSeq) newFree (s->rSeq, s->len, char) ;
      newFree (s->name, strlen(s->name)+1, char) ;
    }
  arrayDestroy (a) ;
}

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: satmatch [options]* <seqfile>\n"
    "  Clusters candidate repeat unit sequences by alignment, in both orientations, of each sequence\n"
    "  against the shorter one of the pair as a circle, which can wrap round any number of times, and\n"
    "  picks a representative for each cluster. By default sequences are taken from the shortest, and\n"
    "  each joins the cluster of the representative so far that it matches best, or starts a new one;\n"
    "  then clusters whose representatives score above -minScore against each other, aligned globally with\n"
    "  satscan's scoring, are merged by single linkage.\n"
    "  Options:\n"
    "    -fa <file_name>         write the representative of each cluster to this fasta file\n"
    "    -1 <file_name>          write the representative of each cluster of more than one sequence to\n"
    "                            this binary ONEcode .1lib file, with type (Y) Satellite\n"
    "    -sf <file_name>         write the score of each pair of sequences to this tab-separated file,\n"
    "                            only with -allPairs\n"
    "    -ti <n>                 transition penalty [1]\n"
    "    -tv <n>                 transversion penalty [2]\n"
    "    -gap <n>                gap penalty per position [3]\n"
    "    -minlen <n>             only use sequences at least this long [0]\n"
    "    -random <n> <len>       use n random sequences of length len instead of <seqfile>\n"
    "    -allPairs               cluster by single linkage over all pairs of sequences\n"
    "    -minScore <n>           link representatives whose satscan score is above this [0]\n"
    "    -T <n>                  number of threads [8]\n" ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit(1) ; }

  FILE *scoreFile = 0, *fastaFile = 0 ;
  char *libFileName = 0, *seqFileName = 0 ;
  ScoreParams scoreParams = {1,2,3} ;
  int nRandom = 0, lenRandom = 0 ;
  int minSeqLen = 0 ;
  int nThread = 8 ;
  int minScore = 0 ;
  bool isAllPairs = false ;

  /* Process command line options */
  while (argc > 1 && **argv == '-')
    { if (!strcmp (*argv, "-sf") && argc > 1)
        { if (!(scoreFile = fopen (*++argv, "w")))
            { fprintf (stderr, "failed to open output file %s\n", *argv) ;
              exit (1) ;
            }
          argc-- ;
        }
      else if (!strcmp (*argv, "-ti") && argc > 1)
        { scoreParams.transition = atoi (*++argv) ;
          if (scoreParams.transition <= 0) die ("bad transition score %s\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-tv") && argc > 1)
        { scoreParams.transversion = atoi (*++argv) ;
          if (scoreParams.transversion <= 0) die ("bad transversion score %s\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-gap") && argc > 1)
        { scoreParams.gap = atoi (*++argv) ;
          if (scoreParams.gap <= 0) die ("bad gap score %s\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-minlen") && argc > 1)
        { minSeqLen = atoi (*++argv) ;
          if (minSeqLen < 0) die ("bad minlen %s - must be non-negative\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-minScore") && argc > 1)
        { minScore = atoi (*++argv) ; argc-- ; }
      else if (!strcmp (*argv, "-allPairs"))
	isAllPairs = true ;
      else if (!strcmp (*argv, "-T") && argc > 1)
        { nThread = atoi (*++argv) ;
          if (nThread < 1) die ("bad number of threads %s - must be at least 1\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp(*argv, "-random") && argc > 2)
	{ nRandom = atoi (*++argv) ;
	  lenRandom = atoi (*++argv) ;
	  argc -= 2 ;
	}
      else if (!strcmp (*argv, "-1") && argc > 2)
	{ libFileName = *++argv ; argc-- ; }
      else if (!strcmp (*argv, "-fa") && argc > 2)
	{ fastaFile = fopen (*++argv, "w") ;
	  if (!fastaFile) die ("failed to open fasta file %s to write", *--argv) ;
	  argc-- ;
	}
      else die ("unknown option: %s\n%s", *argv, usage) ;
      argc-- ; argv++ ;
    }
  
  if (!(argc == 1 || (nRandom > 0 && argc == 0))) die (usage) ;
  if (scoreFile && !isAllPairs) die ("-sf needs -allPairs\n%s", usage) ;
   
  Array aSeq ; 
  int totSeq = 0, nShort = 0 ;
  if (nRandom > 0)
    { aSeq = randomSeqs (nRandom, lenRandom) ;
      totSeq = nRandom * lenRandom ;
    }
  else
    { seqFileName = *argv ;
      SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
      if (!sio) die ("failed to open sequence file %s", *argv) ;
      aSeq = arrayCreate (256, Seq) ;
      while (seqIOread (sio))
	{ if (sio->seqLen < minSeqLen) { ++nShort ; continue ; }
	  Seq *s = arrayp (aSeq, arrayMax(aSeq), Seq) ;
	  s->len = sio->seqLen ;
	  s->name = strdup (sqioId(sio)) ;
	  s->seq = new (s->len, char) ; memcpy (s->seq, sqioSeq(sio), s->len) ;
	  s->rSeq = seqRevComp (s->seq, s->len) ;
	  s->level = satSchemeLevel (s->seq, s->len) ;
	  totSeq += s->len ;
	}
      seqIOclose (sio) ;
    }
  fprintf (stderr, "%d sequences, total length %d", (int)arrayMax(aSeq), totSeq) ;
  if (nShort) fprintf (stderr, ", after skipping %d shorter than %d", nShort, minSeqLen) ;
  fputc ('\n', stderr) ;
  int nSeq = arrayMax(aSeq) ;

  // estimate significance threshold
  double threshFac ;
  { static int NSEQ = 40, LEN = 100 ;
    Array a = randomSeqs (NSEQ, LEN) ;
    int sum = 0 ;
    for (int i = 0 ; i < NSEQ ; ++i)
      for (int j = i+1 ; j < NSEQ ; ++j)
	sum += needlemanWunsch (arrp(a,i,Seq), arrp(a,j,Seq), &scoreParams, NO_THRESHOLD, 0) ;
    seqArrayDestroy (a) ;
    threshFac = 0.75 * 0.01 * sum * 2.0 / (NSEQ*(NSEQ-1)) ; // a hack, but seems good
    fprintf (stderr, "threshold %.3f\n", threshFac) ;
  }
  
  int *cluster = new0 (nSeq, int) ; // 1-based cluster number of each sequence
  int nCluster = 0 ;
  int *scoreij = 0 ;	  // with -allPairs, the scores of all pairs i < j
  if (scoreFile) fprintf (scoreFile, "#seq1\tlen1\tseq2\tlen2\tscore\torientation\toffset\n") ;
  if (isAllPairs)
    { // score all the pairs, in parallel
      scoreij = new0 ((size_t)nSeq*nSeq, int) ;
      PairJob pj = { aSeq, &scoreParams, nSeq, scoreij, scoreFile ? new0 ((size_t)nSeq*nSeq, int) : 0,
		     scoreFile ? 0 : threshFac, 0 } ; // with -sf the scores of all pairs are wanted
      pthread_mutex_init (&pj.lock, 0) ;
      runThreads (pairThread, &pj, nThread) ;
      pthread_mutex_destroy (&pj.lock) ;
      fprintf (stderr, "scoring all pairs\t") ; timeUpdate (stderr) ;

      // now go through the pairs in order, linking them into clusters by single linkage
      for (int i = 0 ; i < nSeq ; i++)
	for (int j = i+1 ; j < nSeq ; j++)
	  { Seq *si = arrp(aSeq,i,Seq) ;
	    Seq *sj = arrp(aSeq,j,Seq) ;
	    int score = scoreij[(size_t)i*nSeq+j] ;
	    if (scoreFile)
	      { int offset = pj.offsetij[(size_t)i*nSeq+j] ;
		fprintf (scoreFile, "%s\t%d\t%s\t%d\t%d\t%c\t%d\n", si->name, si->len, sj->name, sj->len,
			 score, offset >= 0 ? '+' : '-', offset >= 0 ? offset : -offset-1) ;
	      }
	    int threshold = threshFac * (si->len < sj->len ? sj->len : si->len) - 10 ;
	    if (score <= threshold) // link it
	      { if (cluster[i] && cluster[j])
		  { if (cluster[i] != cluster[j])
		      { int from, to ; to = min (cluster[i], cluster[j]) ; from = cluster[i]+cluster[j] - to ;
			printf ("  merge mark %d to %d via %d (%d) to %d (%d) score %d thresh %d\n",
				from, to, i, si->len, j, sj->len, score, threshold) ;
			for (int k = 0 ; k < nSeq ; ++k) if (cluster[k] == from) cluster[k] = to ;
		      }
		  }
		else if (!cluster[i] && !cluster[j])
		  { cluster[i] = cluster[j] = ++nCluster ;
		    printf ("  new mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			    nCluster, i, si->len, j, sj->len, score, threshold) ;
		  }
		else if (cluster[i] && !cluster[j])
		  { cluster[j] = cluster[i] ;
		    printf ("  link mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			    cluster[i], i, si->len, j, sj->len, score, threshold) ;
		  }
		else if (!cluster[i] && cluster[j])
		  { cluster[i] = cluster[j] ;
		    printf ("  link mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			    cluster[j], j, sj->len, i, si->len, score, threshold) ;
		  }
	      }
	  }
      if (pj.offsetij) newFree (pj.offsetij, (size_t)nSeq*nSeq, int) ;
    }
  else
    nCluster = clusterGreedy (aSeq, &scoreParams, threshFac, nThread, cluster) ;

  int *minLen = new (nCluster+1, int), *rep = new (nCluster+1, int) ;
  int *repScore = new0 (nSeq, int) ; // score of each sequence to the representative of its cluster
  chooseRepresentatives (aSeq, &scoreParams, nThread, cluster, nCluster, 0, scoreij, scoreFile != 0,
			 threshFac, minLen, rep, repScore) ;
  if (!isAllPairs) // link the clusters via their representatives, and choose again for merged ones
    { bool *isMerged = new0 (nCluster+1, bool) ;
      if (linkRepresentatives (aSeq, &scoreParams, minScore, nThread, cluster, nCluster, rep, isMerged))
	chooseRepresentatives (aSeq, &scoreParams, nThread, cluster, nCluster, isMerged, 0, false,
			       threshFac, minLen, rep, repScore) ;
      newFree (isMerged, nCluster+1, bool) ;
    }
  Array *member = new0 (nCluster+1, Array) ; // the members of each cluster, in input order
  for (int i = 0 ; i < nSeq ; ++i)
    if (cluster[i])
      { if (!member[cluster[i]]) member[cluster[i]] = arrayCreate (8, int) ;
	array(member[cluster[i]], arrayMax(member[cluster[i]]), int) = i ;
      }

  int nUnsafe = 0 ; // representatives safe for no scheme, left out of -fa and -1
  OneFile *ofLib = 0 ;
  OneSchema *schema = 0 ;
  char *libSeq = 0 ; int libSeqMax = 0 ;
  if (libFileName)
    { schema = oneSchemaCreateFromText (repeatSchemaText) ;
      if (!(ofLib = oneFileOpenWriteNew (libFileName, schema, "lib", true, 1)))
	die ("failed to open .1lib file %s to write", libFileName) ;
      oneAddProvenance (ofLib, "satmatch", VERSION, getCommandLine()) ;
      if (!nRandom) oneAddReference (ofLib, seqFileName, 1) ;
    }
  for (int k = 1 ; k <= nCluster ; ++k)
    { if (!member[k]) continue ;
      Seq *s = arrp(aSeq, rep[k], Seq) ;
      char id[64] ; sprintf (id, "sat-%d-%d-%d", k, s->len, n50Copies (aSeq, member[k])) ;
      bool isSafe = s->level < MAX_SCHEME ; // else satscan would wall its own arrays off, so it is dropped
      if (!isSafe) { fprintf (stderr, "%s%s", nUnsafe ? " " : "dropped as unsafe for every scheme: ", id) ; ++nUnsafe ; }
      if (ofLib && isSafe && arrayMax(member[k]) > 1)
	{ if (s->len > libSeqMax)
	    { if (libSeq) newFree (libSeq, libSeqMax, char) ;
	      libSeqMax = s->len ; libSeq = new (libSeqMax, char) ;
	    }
	  for (int i = 0 ; i < s->len ; ++i) libSeq[i] = index2char[s->seq[i]] ;
	  oneWriteLine (ofLib, 'S', s->len, libSeq) ;
	  oneWriteLine (ofLib, 'I', strlen(id), id) ;
	  oneWriteLine (ofLib, 'Y', 9, "Satellite") ;
	}
      printf ("cluster %d count %d minLen %d repLen %d name %s ", k, (int)arrayMax(member[k]), minLen[k], s->len, id) ;
      for (int i = 0 ; i < s->len ; ++i) putchar (index2char[s->seq[i]]) ;
      putchar ('\n') ;
      if (fastaFile && isSafe)
	{ fprintf (fastaFile, ">%s\n", id) ;
	  for (int i = 0 ; i < s->len ; ++i) fputc (index2char[s->seq[i]], fastaFile) ;
	  fputc ('\n', fastaFile) ;
	}
      for (int a = 0 ; a < arrayMax(member[k]) ; ++a)
	{ Seq *s = arrp(aSeq, arr(member[k],a,int), Seq) ;
	  printf ("  %s\t%d\t%d\n", s->name, s->len, repScore[arr(member[k],a,int)]) ;
	}
    }

  for (int k = 1 ; k <= nCluster ; ++k) if (member[k]) arrayDestroy (member[k]) ;
  newFree (member, nCluster+1, Array) ;
  newFree (minLen, nCluster+1, int) ; newFree (rep, nCluster+1, int) ;
  newFree (repScore, nSeq, int) ;
  newFree (cluster, nSeq, int) ;
  if (scoreij) newFree (scoreij, (size_t)nSeq*nSeq, int) ;
  if (nUnsafe) fputc ('\n', stderr) ;
  if (ofLib) { oneFileClose (ofLib) ; oneSchemaDestroy (schema) ; } // close writes the footer
  if (libSeq) newFree (libSeq, libSeqMax, char) ;
  if (fastaFile) fclose (fastaFile) ;
  if (scoreFile) fclose (scoreFile) ;
  seqArrayDestroy (aSeq) ;
  timeTotal (stderr) ;
  return 0 ;
}

/* Needleman-Wunsch alternative implementation without conditional tests, for fast pipelining
   
  Fixed 1 for transition, 2 for transversion, 3 for gap

  Standard update from row z0[j=0..n-1] to z1[j=0..n-1] is
    z1[j+1] = min (z1[j]+3,z0[j]+v[j+1],z0[j+1]+3) = min (z0[j]+v[j+1], 3+min(z1[j], z0[j+1]))
  Here v[j] is the score for matching x[j] to the character for the current row. We can prebuild
  four vectors v0[], v1[], v2[], v[3] for these depending on the character.

  Alternative is to update a vector of values d1[j] = (3+diff to above)<<3 | (3+diff to left)
  NB diff to left (or above) must be in range -3..3, so (3+diff) is in range 0..6.
  Update is d1[j+1] = f(v[j+1], d1[j] & 0x38, d0[j] & 0x03) = nwMap(w[j] | d1[j]&0x38 | d0[j]|0x03)
  where w[j] = v[j+1]<<6. So we can store nwMap as a table lookup in a 256-dimension U8 table, and
  then progress with DP by updating d1[] without branching.
*/

U8 *nwMapCreate (void)
{
  U8 *map = new (256, U8) ;

  return map ;
}


/***************************** end of file ******************************/
