
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

typedef struct MatchStruct {
  U64  start ;	   // position in current sequence
  I64  bestScore ; // best score so far at end of a query unit
  U64  bestEnd ;   // end position of bestScore
  I64 *score ;     // length q->len - must be positive
  I64 *gap ;       // length q->len
  int  query ;     // index in aQuery
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
  Match *match ;  // index of first active match in aMatch
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

static void printSeq (char *s)
{ putchar (' ') ; putchar (' ') ;
  for (int i = 0 ; i < 32 ; ++i) putchar ("acgt"[s[i]]) ;
  putchar ('\n') ;
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
    "    -smin <n>          minimum score to report [50]\n" ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit(1) ; }

  FILE *outFile = stdout ;
  ScoreParams sp = {1,-1,-2,-7,-1} ;
  U64 scoreMin = 50 ;

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
  U64 nSeq = 0, totSeq = 0, totStart = 0, totMatch = 0, totMatchLen = 0 ;
  Array aMatch = 0 ;
  while (seqIOread (sio))
    { ++nSeq ; totSeq += sio->seqLen ;
      if (sio->seqLen < 32) continue ;
      U64 nStart = 0, nMatch = 0, matchLen = 0, nActive = 0 ;
      // first make array of starts
      fprintf (stderr,"processing sequence %s length %llu", sqioId(sio), sio->seqLen) ;
      aMatch = arrayReCreate(aMatch, 4096, Match) ; arrayMax(aMatch) = 1 ; // burn first value
      char *s = sqioSeq(sio) + 32 ;
      for (U64 i = 32 ; i < sio->seqLen ; ++i, ++s)
	{ int *bs = baseScore[*s] ;
	  for (int j = 0 ; j < arrayMax(aQuery) ; ++j)
	    { Query *q  = arrp(aQuery, j, Query) ;
	      char  *qs = q->seq ;
	      int newStart = 0 ; for (int k = 0 ; k < 32 ; ++k) newStart += q->startScore[k][s[k-32]] ;
	      int newEnd = 0 ;
	      if (i < sio->seqLen - 32) for (int k = 0 ; k < 32 ; ++k) newEnd += q->endScore[k][s[k]] ;
	      Match *lastMatch = 0, *next ;
	      for (Match *match = q->match ; match ; match = next)
		{
		  next = match->next ; // save before match can be freed below
			  // first rotate
		  I64 *ms = match->score ;
		  I64 msn = ms[q->len-1] ; memmove (ms+1, ms, (q->len-1)*sizeof(I64)) ; ms[0] = msn ;
		  // then add diagonal match scores
		  for (int k = 0 ; k < q->len ; ++k, ++ms) if (*ms) *ms += bs[qs[k]] ;
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
		    { if (match->bestScore > (I64)scoreMin)
			{ fprintf (outFile, "%s\t%llu\t%llu\t%s\t%llu\n", sqioId(sio),
				   match->start, match->bestEnd, q->name, match->bestScore) ;
			  ++nMatch ; matchLen += match->bestEnd - match->start ;
			}
		      if (lastMatch) lastMatch->next = next ; else q->match = next ;
		      newFree (match->score, q->len, I64) ;
		      newFree (match->gap, q->len, I64) ;
		      newFree (match, 1, Match) ;
		      --nActive ;
			      continue ; // match is freed: skip the code below and don't set lastMatch
		    }
		  else
		    {
		      // printf ("  %llu i %llu newEnd %d :", match->start, i, newEnd) ;
		      // for (int k = 0 ; k < q->len ; ++k) printf (" %lld", match->score[k]) ;
		      // putchar ('\n') ;
		      // update bestScore if necessary
		      I64 mEnd = match->score[(32*(q->len-1)) % q->len] ;
		      if ((newEnd > 0) && mEnd && (mEnd + newEnd > match->bestScore))
			{ match->bestScore = mEnd + newEnd ; match->bestEnd = i + 32 ;
			  // printf ("** new best i %llu (32*(q->len-1))%%q->len %d mEnd %lld newEnd %d\n",
			  //  	  i, (32*(q->len-1))%q->len, mEnd, newEnd) ;
			}
		    }
		    
		  if (match->score[32 % q->len] >= newStart) newStart = 0 ; // don't restart a new one

		  lastMatch = match ;
		}
	      
	      if (newStart > 0)
		{ ++nStart ; ++nActive ;
		  Match *match = new0 (1, Match) ;
		  match->query = j ; match->start = i-32 ;
		  match->next = q->match ; q->match = match ; // link into list
		  I64 *ms = match->score = new0 (q->len, I64) ;
		  match->gap = new0 (q->len, I64) ;
		  int k = 32 % q->len ;
		  ms[k] = newStart ;
		  // if (match->next)
		  // printf ("RESTART newStart %lld start %llu last->score %lld last->start %llu\n",
		  // match->score[k], match->start, match->next->score[k], match->next->start) ;
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

	  if (!(i%1000000)) { fputc ('.',stderr) ; fflush (stderr) ; }
	}
      for (int j = 0 ; j < arrayMax(aQuery) ; ++j) // flush matches still active at end of sequence
	{ Query *q = arrp(aQuery, j, Query) ;
	  Match *match = q->match, *next ;
	  for ( ; match ; match = next)
	    { next = match->next ;
	      if (match->bestScore > (I64)scoreMin)
		{ fprintf (outFile, "%s\t%llu\t%llu\t%s\t%llu\n", sqioId(sio),
				   match->start, match->bestEnd, q->name, match->bestScore) ;
		  ++nMatch ; matchLen += match->bestEnd - match->start ;
		}
	      newFree (match->score, q->len, I64) ;
	      newFree (match->gap, q->len, I64) ;
	      newFree (match, 1, Match) ;
	    }
	  q->match = 0 ;
	}
      fprintf (stderr," %llu starts, %llu matches, total match length %llu\n", nStart, nMatch, matchLen) ;
      totStart += nStart ; totMatch += nMatch ; totMatchLen += matchLen ;
    }
  seqIOclose (sio) ;
  fprintf (stderr, "processed %llu target sequences, total length %llu bp\n", nSeq, totSeq) ;
  fprintf (stderr, "found %llu starts, %llu matches, total match length %llu bp\n",
	   totStart, totMatch, totMatchLen) ;
  timeTotal(stderr) ;
  
  queryDestroy (aQuery) ;  
}
