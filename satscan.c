
/*  File: satscan.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description: statellite scanning package
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 12 17:00 2026 (rd109)
 * Created: Sun Aug 16 14:02:47 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "seqio.h"
#include "ONElib.h"
#include <string.h>
#include <time.h>
#include <pthread.h>

typedef struct MatchStruct {
  U64  start ;	   // position in current sequence
  I64  bestScore ; // best score so far at end of a query unit
  U64  bestEnd ;   // end position of bestScore
  I64 *score ;     // length q->len - must be positive
  I64 *gap ;       // length q->len
  int  query ;     // index in aQuery
  int  target ;    // index of the target sequence, in file order
  struct MatchStruct  *next ;      // index of nextMatch in aMatch
  bool isGap ;     // if we need to update gap
} Match ;

typedef struct {
  int   len ;
  int   isRC ; 
  char *name ;
  char *seq ;
  char  startScore[32][4] ;
  char  endScore[32][4] ;
} Query ;

typedef struct {
  int match ;
  int transition ;      /* A<->G, C<->T */
  int transversion ;    /* A<->C, A<->T, G<->C, G<->T */
  int gapOpen ;         /* gap penalty per position */
  int gapExtend ;       /* gap extend per position */
} ScoreParams ;

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

static int baseScore[4][4] ;

static void makeStartScore (Query *q)
{
  for (int i = 0 ; i < 32 ; ++i)
    for (int j = 0 ; j < 4 ; ++j)
      q->startScore[i][j] = baseScore[q->seq[i%q->len]][j] ;
}

static void makeEndScore (Query *q)
{
  for (int i = 0 ; i < 32 ; ++i)
    for (int j = 0 ; j < 4 ; ++j)
      q->endScore[i][j] = baseScore[q->seq[(32*q->len+i-32)%q->len]][j] ;
}

typedef struct {		// read-only once scanning starts, so can be shared between threads
  Array       aQuery ;
  ScoreParams sp ;
  U64         scoreMin ;
} ScanParams ;

typedef struct {		// everything that is private to one scanning thread
  Match **active ;		// [nQuery] heads of the lists of active matches for each query
  Array   aMatch ;		// Match* to report: score and gap are freed, so these are small
  U64     totStart, totMatch, totMatchLen ; // summed over all the targets this worker did
} Worker ;

typedef struct {
  int   target ;		// index of this target sequence, in file order
  char *name ;
  char *seq ;			// not owned by the job
  U64   len ;
  bool  isOwned ;		// true if seq is a malloc'd copy that the worker must free
} Job ;

typedef struct {
  U64    nStart, nMatch, matchLen ;
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
  return w ;
}

static void workerDestroy (Worker *w, ScanParams *scan)
{ newFree (w->active, arrayMax(scan->aQuery), Match*) ;
  arrayDestroy (w->aMatch) ;
  newFree (w, 1, Worker) ;
}

// match is done and unlinked: keep it for reporting if it scores well enough, else free it
static void matchFinish (Worker *w, ScanParams *scan, Match *match, Query *q, ScanStats *st)
{
  newFree (match->score, q->len, I64) ; match->score = 0 ;
  newFree (match->gap, q->len, I64) ; match->gap = 0 ;
  if (match->bestScore > (I64)scan->scoreMin)
    { match->next = 0 ;
      array(w->aMatch, arrayMax(w->aMatch), Match*) = match ;
      ++st->nMatch ; st->matchLen += match->bestEnd - match->start ;
    }
  else
    newFree (match, 1, Match) ;
}

static void scanTarget (ScanParams *scan, Worker *w, Job *job, ScanStats *st)
{
  Array aQuery = scan->aQuery ;
  const ScoreParams sp = scan->sp ;
  int nQuery = arrayMax(aQuery) ;
  memset (st, 0, sizeof(ScanStats)) ;
  double t0 = wallTime () ;

  char *s = job->seq + 32 ;
  for (U64 i = 32 ; i < job->len ; ++i, ++s)
    { int *bs = baseScore[(int)*s] ;
      for (int j = 0 ; j < nQuery ; ++j)
	{ Query *q  = arrp(aQuery, j, Query) ;
	  char  *qs = q->seq ;
	  int newStart = 0 ; for (int k = 0 ; k < 32 ; ++k) newStart += q->startScore[k][(int)s[k-32]] ;
	  int newEnd = 0 ;
	  if (i < job->len - 32) for (int k = 0 ; k < 32 ; ++k) newEnd += q->endScore[k][(int)s[k]] ;
	  Match *lastMatch = 0, *next ;
	  for (Match *match = w->active[j] ; match ; match = next)
	    {
	      next = match->next ; // save before match can be freed below
	      // first rotate
	      I64 *ms = match->score ;
	      I64 msn = ms[q->len-1] ; memmove (ms+1, ms, (q->len-1)*sizeof(I64)) ; ms[0] = msn ;
	      // then add diagonal match scores
	      for (int k = 0 ; k < q->len ; ++k, ++ms) if (*ms) *ms += bs[(int)qs[k]] ;
	      // next right gaps
	      ms = match->score ;
	      I64 d = ms[q->len-1] + sp.gapOpen ;
	      for (int k = 0 ; k < q->len ; ++k, ++ms)
		if (d > *ms) { *ms = d ; d += sp.gapExtend ; }
		else if (*ms + sp.gapOpen > d) d = *ms + sp.gapOpen ;
		else d += sp.gapExtend ;
	      ms = match->score ; // extend wrap-around as far as needed
	      for (int k = 0 ; d > *ms ; ++k)
		{ *ms++ = d ; d += sp.gapExtend ; }
	      // next if necessary do the down gaps
	      if (match->isGap)
		{ ms = match->score ;
		  I64 *mg = match->gap ;
		  match->isGap = false ;
		  for (int k = 0 ; k < q->len ; ++k, ++ms, ++mg)
		    { if (*mg > *ms) { *ms = *mg ; *mg += sp.gapExtend ; }
		      else if (*ms + sp.gapOpen > *mg) *mg = *ms + sp.gapOpen ;
		      else *mg += sp.gapExtend ;
		      if (*mg > 0) match->isGap = true ;
		    }
		}
	      // finally, check if done, and if so then process match
	      ms = match->score ;
	      bool isDone = true ;
	      for (int k = 0 ; k < q->len ; ++k, ++ms)
		if (*ms > 0) isDone = false ;
		else *ms = 0 ;
	      if (isDone)
		{ if (lastMatch) lastMatch->next = next ; else w->active[j] = next ;
		  matchFinish (w, scan, match, q, st) ;
		  continue ; // match is gone: skip the code below and don't set lastMatch
		}
	      else
		{
		  // update bestScore if necessary
		  I64 mEnd = match->score[(32*(q->len-1)) % q->len] ;
		  if ((newEnd > 0) && mEnd && (mEnd + newEnd > match->bestScore))
		    { match->bestScore = mEnd + newEnd ; match->bestEnd = i + 32 ; }
		}

	      if (match->score[32 % q->len] >= newStart) newStart = 0 ; // don't restart a new one

	      lastMatch = match ;
	    }

	  if (newStart > 0)
	    { ++st->nStart ;
	      Match *match = new0 (1, Match) ;
	      match->query = j ; match->target = job->target ; match->start = i-32 ;
	      match->next = w->active[j] ; w->active[j] = match ; // link into list
	      I64 *ms = match->score = new0 (q->len, I64) ;
	      match->gap = new0 (q->len, I64) ;
	      int k = 32 % q->len ;
	      ms[k] = newStart ;
	      int d = newStart + sp.gapOpen ;
	      if (d > 0)
		{ match->isGap = true ;
		  match->gap[k] = d ;
		  while (d > 0)
		    { if (++k == q->len) k = 0 ;
		      if (ms[k] >= d) break ;
		      ms[k] = d ;
		      d += sp.gapExtend ;
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
  w->totStart += st.nStart ; w->totMatch += st.nMatch ; w->totMatchLen += st.matchLen ;
  if (job->isOwned) free (job->seq) ;
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

// sorts by target in file order, then start position, then query
static int matchOrder (const void *a, const void *b)
{
  const Match *x = *(Match* const*)a, *y = *(Match* const*)b ;
  if (x->target != y->target) return x->target < y->target ? -1 : 1 ;
  if (x->start != y->start) return x->start < y->start ? -1 : 1 ;
  if (x->query != y->query) return x->query < y->query ? -1 : 1 ;
  return 0 ;
}

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: satscan <options>* <sat_seqfile> <genome_seqfile>\n"
    "  options are:\n"
    "    -o <file_name>     output file [standard output]\n"
    "    -m <n>             match score [1]\n"
    "    -ti <n>            transition mismatch score [-1]\n"
    "    -tv <n>            transversion mismatch score [-2]\n"
    "    -gapOpen <n>       gap open score [-7]\n"
    "    -gapExtend <n>     gap extend score [-1]\n"
    "    -smin <n>          minimum score to report [50]\n"
    "    -T <n>             number of scanning threads [8]\n" ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit(1) ; }

  FILE *outFile = stdout ;
  ScoreParams sp = {1,-1,-2,-7,-1} ;
  U64 scoreMin = 50 ;
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
      else if (!strcmp (*argv, "-m") && argc > 1)
        { sp.match = atoi (*++argv) ;
          if (sp.match <= 0) die ("bad match score %s - must be +ve\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-ti") && argc > 1)
        { sp.transition= atoi (*++argv) ;
          if (sp.transition >= 0) die ("bad transition score %s - must be -ve\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-tv") && argc > 1)
        { sp.transversion = atoi (*++argv) ;
          if (sp.transversion >= 0) die ("bad transversion score %s - must be -ve\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-gapOpen") && argc > 1)
        { sp.gapOpen = atoi (*++argv) ;
          if (sp.gapOpen >= 0) die ("bad gapOpen score %s\n - must be -ve%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-T") && argc > 1)
        { nThread = atoi (*++argv) ;
          if (nThread < 1) die ("bad number of threads %s - must be at least 1\n%s", *argv, usage) ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-smin") && argc > 1)
        { long long v = atoll (*++argv) ;
          if (v < 0) die ("bad smin score %s - must be non-negative\n%s", *argv, usage) ;
          scoreMin = v ;
          argc-- ;
        }
      else if (!strcmp (*argv, "-gapExtend") && argc > 1)
        { sp.gapExtend = atoi (*++argv) ;
          if (sp.gapExtend >= 0) die ("bad gapExtend score %s\n - must be -ve%s", *argv, usage) ;
          argc-- ;
        }

      else die ("unknown option: %s\n%s", *argv, usage) ;
      argc-- ; argv++ ;
    }
  if (argc != 2) die (usage) ;

  for (int i = 4 ; i-- ;) for (int j = 4 ; j-- ;) baseScore[i][j] = sp.transversion ;
  baseScore[0][0] = baseScore[1][1] = baseScore[2][2] = baseScore[3][3] = sp.match ;
  baseScore[0][2] = baseScore[2][0] = baseScore[1][3] = baseScore[3][1] = sp.transition ;
 
  SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
  if (!sio) die ("failed to open query sequence file %s", *argv) ;
  Array aQuery = arrayCreate (256, Query) ;
  U64 totSat = 0 ;
  while (seqIOread (sio))
    { Query *q = arrayp (aQuery, arrayMax(aQuery), Query) ;
      q->len = sio->seqLen ; 
      q->name = strdup (sqioId(sio)) ;
      q->seq = new (q->len, char) ; memcpy (q->seq, sqioSeq(sio), q->len) ;
      makeStartScore (q) ; makeEndScore (q) ;
      Query *qr = arrayp (aQuery, arrayMax(aQuery), Query) ;
      q = arrp(aQuery, arrayMax(aQuery)-2, Query) ; // arrayp() may have reallocated
      qr->len = q->len ; qr->name = q->name ; qr->isRC = true ;
      qr->seq = seqRevComp (q->seq, q->len) ;
      makeStartScore (qr) ; makeEndScore (qr) ;
      totSat += q->len ;
    }
  seqIOclose (sio) ;
  --argc ; ++argv ;
  fprintf (stderr, "read in %llu query sequences, total length %llu bp\n",
	   arrayMax(aQuery)/2, totSat) ;
 
  sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
  if (!sio) die ("failed to open target sequence file %s", *argv) ;
  ScanParams scan = { aQuery, sp, scoreMin } ;
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
      job.name = strdup (sqioId(sio)) ;
      array(aTargetName, job.target, char*) = job.name ;
      job.len = sio->seqLen ;
      if (nThread == 1)
	{ job.seq = sqioSeq(sio) ; job.isOwned = false ;
	  workerRun (&scan, worker[0], &job) ;
	}
      else
	{ job.seq = malloc (job.len) ;	// copy, because the next seqIOread overwrites the buffer
	  if (!job.seq) die ("failed to allocate %llu bytes for sequence %s", job.len, job.name) ;
	  memcpy (job.seq, sqioSeq(sio), job.len) ; job.isOwned = true ;
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
  U64 totStart = 0, totMatch = 0, totMatchLen = 0 ;
  Array aMatch = worker[0]->aMatch ;
  for (int t = 0 ; t < nThread ; ++t)
    { Worker *w = worker[t] ;
      totStart += w->totStart ; totMatch += w->totMatch ; totMatchLen += w->totMatchLen ;
      if (t)
	for (U64 i = 0 ; i < arrayMax(w->aMatch) ; ++i)
	  array(aMatch, arrayMax(aMatch), Match*) = arr(w->aMatch, i, Match*) ;
    }
  arraySort (aMatch, matchOrder) ;
  for (U64 i = 0 ; i < arrayMax(aMatch) ; ++i)
    { Match *m = arr(aMatch, i, Match*) ;
      fprintf (outFile, "%s\t%llu\t%llu\t%s\t%llu\n", arr(aTargetName, m->target, char*),
	       m->start, m->bestEnd, arrp(aQuery, m->query, Query)->name, m->bestScore) ;
      newFree (m, 1, Match) ;
    }

  fprintf (stderr, "processed %llu target sequences, total length %llu bp\n", nSeq, totSeq) ;
  fprintf (stderr, "found %llu starts, %llu matches, total match length %llu bp\n",
	   totStart, totMatch, totMatchLen) ;
  timeTotal(stderr) ;

  for (U64 i = 0 ; i < arrayMax(aTargetName) ; ++i) free (arr(aTargetName, i, char*)) ;
  arrayDestroy (aTargetName) ;
  for (int t = 0 ; t < nThread ; ++t) workerDestroy (worker[t], &scan) ;
  newFree (worker, nThread, Worker*) ;
  queryDestroy (aQuery) ;
}
