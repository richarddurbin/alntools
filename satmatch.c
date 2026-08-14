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
 * Last edited: Aug 14 16:08 2026 (rd109)
 * Created: Tue Aug 11 08:32:53 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "seqio.h"
#include "ONElib.h"
#include <string.h>

typedef struct {
  int   len ;
  char *name ;
  char *seq ;
  char *rSeq ;
} Seq ;

typedef struct {
  int transition ;      /* A<->G, C<->T */
  int transversion ;    /* A<->C, A<->T, G<->C, G<->T */
  int gap ;             /* gap penalty per position */
  bool isLoop ;
} ScoreParams ;

static inline int min (int x, int y) { return x < y ? x : y ; }

int needlemanWunsch (Seq *s1, Seq *s2, ScoreParams *sp, int *jMax)
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

  // Allocate two rows
  int *prev = new0 (len2+1, int) ; // 0 if looping
  int *row = new (len2+1, int) ;
  // Initialize first row with gaps unless looping
  if (!sp->isLoop) for (j = 0 ; j <= len2 ; j++) prev[j] = j * sp->gap ;

  // Fill DP
  //  for (j = 0 ; j < len1 ; ++j) printf ("  %c", s2->seq[j]) ; putchar ('\n') ;
  for (i = 1 ; i <= len1 ; i++)
    { if (sp->isLoop)
	row[0] = min (prev[len2-1]+cost[s1->seq[i-1]][s2->seq[len2-1]], sp->gap+prev[len2]) ;
      else
	row[0] = i * sp->gap ;
      // printf ("%c %2d", index2char[s1->seq[i-1]], row[0]) ;
      for (j = 1 ; j <= len2 ; j++)
        { int d = prev[j-1] + cost[s1->seq[i-1]][s2->seq[j-1]] ;  // diagonal: match/mismatch
          int g = sp->gap + min (prev[j], row[j-1]) ;             // up or left + gap
          row[j] = min (d, g) ;
	  // printf (" %2d", row[j]) ;
        }
      // putchar ('\n') ;
      int *tmp = prev ; prev = row ; row = tmp ; // swap rows
    }
  
  int score ;
  if (sp->isLoop)
    { score = 3 * len1 ;
      for (j = 0 ; j < len2 ; ++j)
	if (prev[j] < score) { score = prev[j] ; if (jMax) *jMax = j ; }
    }
  else score = prev[len2] ;
  
  newFree (prev, len2+1, int) ;
  newFree (row, len2+1, int) ;

  return score ;
}

Array randomSeqs (int n, int len)
{
  Array a = arrayCreate (n, Seq) ;
  for (int i = 0 ; i < n ; ++i)
    { Seq *s = arrayp (a, arrayMax(a), Seq) ;
      s->len = len ;
      s->name = new(32,char) ; sprintf (s->name, "random%d-%d", len, i+1) ;
      s->seq = new (s->len, char) ; for (int j = 0 ; j < len ; ++j) s->seq[j] = rand()&0x3 ;
    }
  return a ;
}

void seqArrayDestroy (Array a)
{
  for (int i = 0 ; i < arrayMax(a) ; ++i)
    { Seq *s = arrp(a,i,Seq) ;
      newFree (s->seq, s->len, char) ;
      newFree (s->name, strlen(s->name)+1, char) ;
    }
  arrayDestroy (a) ;
}

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: satmatch [-sf <scorefile>] [-fa <fastafile>] [-ti <n>] [-tv <n>] [-gap <n>] [-loop] [-random <n> <len>] <seqfile>" ;
  if (!argc) { fprintf (stderr, "%s\n", usage) ; exit(1) ; }

  FILE *scoreFile = 0, *fastaFile = 0 ;
  ScoreParams scoreParams = {1,2,3,false} ;
  int nRandom = 0, lenRandom = 0 ;

  /* Process command line options */
  while (argc > 1 && **argv == '-')
    { if (!strcmp (*argv, "-sf"))
        { if (!(scoreFile = fopen (*++argv, "w")))
            { fprintf (stderr, "failed to open output file %s\n", *argv) ;
              exit (1) ;
            }
          argc-- ; argv++ ;
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
      else if (!strcmp (*argv, "-loop"))
	{ scoreParams.isLoop = true ; }
      else if (!strcmp(*argv, "-random") && argc > 2)
	{ nRandom = atoi (*++argv) ;
	  lenRandom = atoi (*++argv) ;
	  argc -= 2 ;
	}
      else if (!strcmp (*argv, "-fa") && argc > 2)
	{ fastaFile = fopen (*++argv, "w") ;
	  if (!fastaFile) die ("failed to open fasta file %s to write", *--argv) ;
	  argc-- ;
	}
      else die ("unknown option: %s\n%s", *argv, usage) ;
      argc-- ; argv++ ;
    }
  
  if (!(argc == 1 || (nRandom > 0 && argc == 0))) die (usage) ;

  // estimate significance threshold
  double threshFac ;
  { static int NSEQ = 40, LEN = 100 ;
    Array a = randomSeqs (NSEQ, LEN) ;
    int sum = 0 ;
    for (int i = 0 ; i < NSEQ ; ++i)
      for (int j = i+1 ; j < NSEQ ; ++j)
	sum += needlemanWunsch (arrp(a,i,Seq), arrp(a,j,Seq), &scoreParams, 0) ;
    seqArrayDestroy (a) ;
    threshFac = 0.9 * 0.01 * sum * 2.0 / (NSEQ*(NSEQ-1)) ; // a hack, but seems good
    fprintf (stderr, "threshold %.3f\n", threshFac) ;
  }
   
  Array aSeq ; 
  int totSeq = 0 ;
  if (nRandom > 0)
    { aSeq = randomSeqs (nRandom, lenRandom) ;
      totSeq = nRandom * lenRandom ;
    }
  else
    { SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, "r") ;
      if (!sio) die ("failed to open sequence file %s", *argv) ;
      aSeq = arrayCreate (256, Seq) ;
      while (seqIOread (sio))
	{ Seq *s = arrayp (aSeq, arrayMax(aSeq), Seq) ;
	  s->len = sio->seqLen ;
	  s->name = strdup (sqioId(sio)) ;
	  s->seq = new (s->len, char) ; memcpy (s->seq, sqioSeq(sio), s->len) ;
	  s->rSeq = seqRevComp (s->seq, s->len) ;
	  totSeq += s->len ;
	}
      seqIOclose (sio) ;
    }
  
  fprintf (stderr, "%d sequences, total length %d\n", (int)arrayMax(aSeq), totSeq) ;

  int nSeq = arrayMax(aSeq) ;
  int *mark = new0 (nSeq, int) ;
  int *scoreij = new0 (nSeq*nSeq, int) ;
  int nMark = 0 ;
  if (scoreFile) fprintf (scoreFile, "#seq1\tlen1\tseq2\tlen2\tscore\n") ;
  for (int i = 0 ; i < nSeq ; i++)
    for (int j = i+1 ; j < nSeq ; j++)
      { Seq *si = arrp(aSeq,i,Seq) ;
        Seq *sj = arrp(aSeq,j,Seq) ;
	int jMax1 = 0, jMax2 = 0 ;
        int score1 = needlemanWunsch(si, sj, &scoreParams, &jMax1) ;
	bool isRC = false ;
	char *tmp = sj->seq ; sj->seq = sj->rSeq ; sj->rSeq = tmp ;
	int score2 = needlemanWunsch(si, sj, &scoreParams, &jMax2) ;
	tmp = sj->seq ; sj->seq = sj->rSeq ; sj->rSeq = tmp ;

	int score = min (score1, score2) ;
	scoreij[i*nSeq+j] = score ;
	if (scoreFile)
	  { if (score == score1)
	      fprintf (scoreFile, "%s\t%d\t%s\t%d\t%d\t+\t%d\n",
		       si->name, si->len, sj->name, sj->len, score1, jMax1) ;
	    else
	      fprintf (scoreFile, "%s\t%d\t%s\t%d\t%d\t-\t%d\n",
		       si->name, si->len, sj->name, sj->len, score2, jMax2) ;
	  }
	
	int threshold = threshFac * (si->len < sj->len ? sj->len : si->len) - 10 ;
	if (score <= threshold) // link it
	  { if (mark[i] && mark[j])
	      { if (mark[i] != mark[j])
		  { int from, to ; to = min (mark[i], mark[j]) ; from = mark[i]+mark[j] - to ;
		    printf ("  merge mark %d to %d via %d (%d) to %d (%d) score %d thresh %d\n",
			    from, to, i, si->len, j, sj->len, score, threshold) ;
		    for (int k = 0 ; k < nSeq ; ++k) if (mark[k] == from) mark[k] = to ;
		  }
	      }
	    else if (!mark[i] && !mark[j])
	      { mark[i] = mark[j] = ++nMark ;
		printf ("  new mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			nMark, i, si->len, j, sj->len, score, threshold) ;
	      }
	    else if (mark[i] && !mark[j])
	      { mark[j] = mark[i] ;
		printf ("  link mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			mark[i], i, si->len, j, sj->len, score, threshold) ;
	      }
	    else if (!mark[i] && mark[j])
	      { mark[i] = mark[j] ;
		printf ("  link mark %d: %d (%d) to %d (%d) score %d thresh %d\n",
			mark[j], j, sj->len, i, si->len, score, threshold) ;
	      }
	  }
      }

  int *sumScore = new0 (nSeq, int) ;
  for (int i = 0 ; i < nSeq ; ++i)
    for (int j = 0 ; j < nSeq ; ++j)
      if (mark[i] && mark[i] == mark[j])
	{ sumScore[i] += scoreij[i*nSeq+j] ;
	  sumScore[j] += scoreij[i*nSeq+j] ;
	}

  for (int k = 1 ; k < nMark ; ++k)
    { int count = 0, minLen = 1<<30, minSum = 1<<30, iMinSum = -1 ;
      for (int i = 0 ; i < nSeq ; ++i)
	if (mark[i] == k)
	  { ++count ;
	    Seq *s = arrp(aSeq,i,Seq) ;
	    if (s->len < minLen) minLen = s->len ;
	    if (sumScore[i] < minSum) { minSum = sumScore[i] ; iMinSum = i ; }
	  }
      if (!count) continue ;
      if (arrp(aSeq,iMinSum,Seq)->len > 1.1*minLen) // need to look again
	{ minSum = 1 << 30 ; iMinSum = -1 ;
	  for (int i = 0 ; i < nSeq ; ++i)
	    if (mark[i] == k)
	      { Seq *s = arrp(aSeq,i,Seq) ;
		if (s->len <= 1.1*minLen && sumScore[i] < minSum)
		  { minSum = sumScore[i] ; iMinSum = i ; }
	      }
	}
      Seq *s = arrp(aSeq, iMinSum, Seq) ;
      printf ("cluster %d count %d minLen %d repLen %d ", k, count, minLen, s->len) ;
      for (int i = 0 ; i < s->len ; ++i) putchar (index2char[s->seq[i]]) ;
      putchar ('\n') ;
      if (fastaFile)
	{ fprintf (fastaFile, ">sat-%d-%d\n", k, s->len) ;
	  for (int i = 0 ; i < s->len ; ++i) fputc (index2char[s->seq[i]], fastaFile) ;
	  fputc ('\n', fastaFile) ;
	}
      for (int i = 0 ; i < nSeq ; ++i)
	if (mark[i] == k)
	  { Seq *s = arrp(aSeq,i,Seq) ;
	    printf ("  %s\t%d\t%d\n", s->name, s->len, scoreij[i*nSeq+iMinSum]) ;
	  }
    }
 
  if (fastaFile) fclose(fastaFile) ;
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
