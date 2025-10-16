/*  File: tanbed.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: make a bed file from a FasTAN .1aln file
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 16 05:10 2025 (rd109)
 * Created: Thu Oct 16 02:48:43 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"

Gdb *readGdb (OneFile *of, int k)
{
  if (!oneGoto (of, 'g', k)) die ("failed to go to GDB %d in %s", k, oneFileName (of)) ;
  oneReadLine (of) ; // swallow the g line

  Gdb *gdb = new0 (1, Gdb) ;
  int N = oneCountUntilNext (of, 'S', 'g') ;
  gdb->seqDict = dictCreate (N) ;
  gdb->seqLen = new0 (N, I64) ;
  gdb->child = new0 (N, int) ;
  // gdb->nCtg = oneCountUntilNext (of, 'C', 'g') ; // 'C' not an object type
  I64 maxCount ; oneStatsContains (of, 'g', 'C', &maxCount, 0) ;
  gdb->ctgLen = new0 (maxCount, I64) ;
  gdb->parent = new0 (maxCount, int) ;
  gdb->offset = new0 (maxCount, I64) ;

  oneReadLine (of) ; // read the first S line, we hope!
  while (of->lineType == 'S')
    { dictAdd (gdb->seqDict, oneString(of), 0) ;
      gdb->child[gdb->nSeq] = gdb->nCtg ;
      I64 end = 0 ;
      while (oneReadLine (of))
	{ bool doneSeq = false ;
	  switch (of->lineType)
	    {
	    case 'S': case 'A': case 'g': doneSeq = true ; break ;
	    case 'G': end += oneInt(of,0) ; break ;
	    case 'C':
	      gdb->parent[gdb->nCtg] = gdb->nSeq ;
	      gdb->ctgLen[gdb->nCtg] = oneInt(of,0) ;
	      gdb->offset[gdb->nCtg] = end ;
	      end += oneInt(of,0) ;
	      ++gdb->nCtg ;
	      break ;
	    case 'M': break ;
	    }
	  if (doneSeq) break ;
	}
      gdb->seqLen[gdb->nSeq] = end ;
      if (gdb->child[gdb->nSeq] == gdb->nCtg) gdb->child[gdb->nSeq] = -1 ; // no contigs
      ++gdb->nSeq ;
    }
  return gdb ;
}

static inline int ctg2seq (Gdb *gdb, int ctg) { return gdb->parent[ctg] ; }

static inline I64 ctg2pos (Gdb *gdb, int ctg, I64 x) { return gdb->offset[ctg] + x ; }

typedef struct {
  int seq ;
  I64 start, end ;
  int unit ;
  int score ;
} BedLine ;

static int bedSort (const void *a, const void *b)
{
  BedLine *ba = (BedLine*)a, *bb = (BedLine*)b ;
  if (ba->seq != bb->seq) return ba->seq - bb->seq ;
  if (ba->start != bb->start) return ba->start - bb->start ;
  return ba->end - bb->end ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (argc != 1) die ("Usage: tanbed <.1aln file>") ;
  OneSchema *schema = oneSchemaCreateFromText (alnSchemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "aln", 1) ;
  if (!of) die ("failed to open .1aln file %s", *argv) ;
  Gdb *gdb = readGdb (of, 1) ;
  if (of->lineType != 'A') die ("unexpected line type %c", of->lineType) ;
  I64 nAlign = 0 ;
  oneStats (of, 'A', &nAlign, 0, 0) ;
  Array ab = arrayCreate (nAlign, BedLine) ;
  I64 totAlign = 0 ;
  while (of->lineType == 'A')
    { BedLine *b = arrayp(ab, arrayMax(ab), BedLine) ;
      int ctg = oneInt(of,0) ;
      if (ctg != oneInt(of,3))
	die ("target mismatch line %lld - not a TAN file?", (long long)of->line) ;
      b->seq = ctg2seq(gdb,ctg) ;
      b->start = ctg2pos(gdb,ctg,oneInt(of,4)) ;
      b->end = ctg2pos(gdb,ctg,oneInt(of,2)) ;
      double len = oneInt(of,2) - oneInt(of,4) ;
      totAlign += oneInt(of,2) - oneInt(of,4) ;
      while (oneReadLine(of) && of->lineType != 'A')
	if (of->lineType == 'D') b->score = (int)(1000*(1.0 - oneInt(of,0)/len)) ;
	else if (of->lineType == 'U') b->unit = oneInt(of,0) ;
    }
  arraySort (ab, bedSort) ;
  int i ;
  for (i = 0 ; i < arrayMax(ab) ; ++i)
    { BedLine *b = arrp(ab, i, BedLine) ;
      printf ("%s\t", dictName(gdb->seqDict, b->seq)) ;
      printf ("%lld\t%lld\t%d\t%d\n", (long long)b->start, (long long)b->end, b->unit, b->score) ;
    }
  I64 totSeq = 0 ; for (i = 0 ; i < gdb->nSeq ; ++i) totSeq += gdb->seqLen[i] ;
  fprintf (stderr, "processed %lld alignments total length %lld from %s length %lld (%.1f %%)\n",
	   nAlign, totAlign, *argv, totSeq, totAlign/(0.01*totSeq)) ;
}
