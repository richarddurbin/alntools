
/*  File: tancons.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: make a bed file from a FasTAN .1aln file
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 10:12 2025 (rd109)
 * Created: Thu Aug 7 02:48:43 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"
#include "seqio.h"

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;
  if (argc != 4) die ("Usage: tancons <unit> <.1aln file> <.1ano> <seqfile>") ;

  int unit = atoi(*argv) ;
  --argc ; ++argv ;

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "aln", 1) ;
  if (!of) die ("failed to open .1aln file %s", *argv) ;
  Gdb *gdb = readGdb (of, 1, stderr) ;
  if (of->lineType != 'A') die ("unexpected line type %c", of->lineType) ;
  --argc ; ++argv ; 

  OneFile *of2 = oneFileOpenRead (*argv, schema, "ano", 1) ;
  if (!of2) die ("failed to open .1ano file %s", *argv) ;
  Gdb *gdb2 = readGdb (of2, 1, stderr) ;
  if (of2->lineType != 'M') die ("unexpected line type %c", of2->lineType) ;
  if (gdb->nSeq != gdb2->nSeq || gdb->nCtg != gdb2->nCtg || gdb->totSeq != gdb2->totSeq || gdb->totCtg != gdb2->totCtg) 
    die ("number of sequences or columns in .1aln and .1ano files differ") ;
  --argc ; ++argv ;

  SeqIO *sio = seqIOopenRead (*argv, dna2indexConv, "r") ;
  if (!sio) die ("failed to open sequence file %s", *argv) ;
  
  I64 nAlign = 0 ;
  I64 uMax = 0, maxLength = 0, maxSeq, maxStart ;
  while (of->lineType == 'A')
    { ++nAlign ;
      int ctg = oneInt(of,0) ;
      if (ctg != oneInt(of,3))
	die ("target mismatch line %lld - not a TAN file?", (long long)of->line) ;
      I64 start = ctg2pos(gdb,ctg,oneInt(of,4)) ;
      I64 len = oneInt(of,2) - oneInt(of,4) ;
      int u = 0 ;
      while (oneReadLine(of) && of->lineType != 'A')
	if (of->lineType == 'U') u = oneInt(of,0) ;
      if (u == unit && len > maxLength) 
        { uMax = nAlign ; maxLength = len ; 
          maxSeq   = ctg2seq(gdb,ctg) ;
          maxStart = start ;
        }
    }
  oneFileClose (of) ;
  char *maxSeqName = dictName(gdb->seqDict, maxSeq) ;
  fprintf (stderr,"longest tandem repeat of unit %d is %lld bases in sequence %s starting at %lld\n", 
          unit, (long long)maxLength, maxSeqName, (long long)maxStart) ;
  
  I64 nParse, nStart = 0, *starts ;
  if (!oneGoto (of2, 'M', 1)) die ("can't locate to first M line in .1ano file") ;
  while (oneReadLine (of2))
    if (of2->lineType == 'M' && oneInt(of2,0) == maxSeq && oneInt(of2,1) == maxStart)
      while (oneReadLine (of2))
        { if (of2->lineType == 'L' && unit != atoi (oneString(of2)))
            die ("unit mismatch %s != %d\n", oneString(of2), unit) ;
          else if (of2->lineType == 'P')
            { nParse = oneLen(of2) ;
              I64 *parse = oneIntList(of2) ;
              starts = new (nParse, I64) ;
              for (int i = 0 ; i < nParse-1 ; ++i)
                if (parse[i+1] - parse[i] == unit) starts[nStart++] = maxStart + parse[i] ;
              break ; // we are done with the longest match, no need to read further
            }
          else if (of2->lineType == 'M') die ("missing P line in .1ano file") ;
        }
  oneFileClose (of2) ;
  fprintf (stderr,"found %lld starts of exact unit size, ", (long long)nStart) ;
 
  // now find the starts in the sequence file and build the consensus
  int i, j ;
  int counts[unit][5], sumCount = 0 ; // A,C,G,T,N
  memset (counts, 0, unit*5*sizeof(int)) ;
  char seq[unit] ;
  while (seqIOread(sio))
    if (!strcmp(sqioId(sio), maxSeqName))
      { if (starts[nStart-1] + unit > sio->seqLen) 
          die ("start %d + unit %d exceeds sequence length %d", (int)starts[nStart-1], unit, (int)sio->seqLen) ; ;
        for (i = 0 ; i < nStart ; ++i) for (j = 0 ; j < unit ; ++j) counts[j][sqioSeq(sio)[starts[i]+j]]++ ;
        for (i = 0 ; i < unit ; ++i)
          { int maxCount = 0, maxBase = 4 ;
            for (j = 0 ; j < 4 ; ++j) if (counts[i][j] > maxCount) { maxCount = counts[i][j] ; maxBase = j ; }
            seq[i] = index2char[maxBase] ;
            sumCount += maxCount ;
          }
        break ;
      }
  seqIOclose (sio) ;
  fprintf (stderr,"with average_match %.1f%%\n", 100.0*sumCount/(nStart*unit)) ; 
  printf (">cons%d\n%s\n", unit, seq) ;
  newFree (starts, nParse, I64) ;
}