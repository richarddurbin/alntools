/*  File: taco.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: TAndem COmpressor
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 12:57 2025 (rd109)
 * Created: Sat Nov 15 23:37:43 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"
#include "seqio.h"

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  char *outFileName = 0 ;
  if (argc >= 2 && !strcmp(*argv, "-o"))
    { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
  
  if (argc != 2)
    { fprintf (stderr, "Usage: taco [-o <outFileName>] <input.1aln> <seqFile>\n"
	       "  taco stands for 'TAndem COmpress (cf hoco for 'HOmopolymer COmpress'\n"
	       "  input.1aln should be created by FasTAN and the names and lengths must match to seqFile\n"
	       "  default outFileName is <seqFile-stem>-taco.1seq;"
	       "  user given outFileName can end in .fa or .fa.gz or .1seq\n") ;
      exit (1) ;
    }
  
  // open the input .1aln file
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *ofIn = oneFileOpenRead (argv[0],schema, "aln", 1) ;
  if (!ofIn) die ("failed to open %s as a .1aln file", argv[0]) ;
  Gdb *gdb = readGdb (ofIn, 1, stderr) ;
  if (ofIn->lineType != 'A') die ("unexpected line type %c", ofIn->lineType) ;

  // and the input sequence
  SeqIO *inIO = seqIOopenRead (argv[1], dna2textConv, false) ;
  if (!inIO) die ("failed to open %s to read as a sequence file", argv[1]) ;

  // and the output file
  if (!outFileName)
    { outFileName = new0 (strlen(argv[1]) + 10, char) ;
      strcpy (outFileName, argv[1]) ;
      char *s = outFileName ; while (*s) ++s ;
      if (s > outFileName+3 && !strcmp(s-3, ".gz")) { s -= 3 ; *s = 0 ; } // remove ".gz"
      while (s > outFileName && *s != '.') --s ;
      if (s == outFileName) s += strlen (outFileName) ; else *s = 0 ; // remove ending from '.'
      strcat (outFileName, "-taco.fa.gz") ;
    }
  SeqIO *outIO = seqIOopenWrite (outFileName, 0, dna2textConv, 0) ;
  if (!outIO) die ("failed to open %s to write as a sequence file", outFileName) ;

  // now read the .1aln file
  I64 nAlign = 0 ;
  oneStats (ofIn, 'A', &nAlign, 0, 0) ;
  Array at = arrayCreate (nAlign, TanLine) ;
  I64 totAlign = 0 ;
  while (ofIn->lineType == 'A')
    { TanLine *t = arrayp(at, arrayMax(at), TanLine) ;
      int ctg = oneInt(ofIn,0) ;
      if (ctg != oneInt(ofIn,3))
	die ("target mismatch line %lld - not a TAN file?", (long long)ofIn->line) ;
      t->seq     = ctg2seq(gdb,ctg) ;
      t->start   = ctg2pos(gdb,ctg,oneInt(ofIn,4)) ;
      t->end     = ctg2pos(gdb,ctg,oneInt(ofIn,2)) ;
      double len = oneInt(ofIn,2) - oneInt(ofIn,4) ;
      totAlign  += oneInt(ofIn,2) - oneInt(ofIn,4) ;
      while (oneReadLine(ofIn) && ofIn->lineType != 'A')
	if (ofIn->lineType == 'D') t->score = (int)(1000*(1.0 - oneInt(ofIn,0)/len)) ;
	else if (ofIn->lineType == 'U') t->unit = oneInt(ofIn,0) ;
    }
  oneFileClose (ofIn) ;
  arraySort (at, tanSort) ;

  // now build a lookup from the sequence name to start position in at
  int i ;
  I32 *seqStart = new0 (gdb->nSeq, I32) ;
  I32 lastSeq = -1 ;
  for (i = 0 ; i < arrayMax(at) ; ++i)
    { I32 seq = arrp(at, i, TanLine)->seq ;
      if (seq != lastSeq) { seqStart[seq] = i ; lastSeq = seq ; }
    } // NB sequences with no matches will have value 0 - that works

  // and make the buffer for writing the tandem-compressed sequence
  I64 maxSeqLen = 0 ;
  for (i = 0 ; i < gdb->nSeq ; ++i)
    if (gdb->seqLen[i] > maxSeqLen)
      maxSeqLen = gdb->seqLen[i] ;
  char *tacoSeq = new (maxSeqLen, char) ;

  // now run through inIO, writing out outIO
  while (seqIOread (inIO))
    { U32 seq ;
      char *id = sqioId(inIO) ;
      if (inIO->descLen) // need to deal with Gene including the description in the ID
	{ id = malloc (strlen(id) + strlen(sqioDesc(inIO)) + 2) ;
	  sprintf (id, "%s %s", sqioId(inIO), sqioDesc(inIO)) ;
	}
      if (!dictFind (gdb->seqDict, id, &seq))
	die ("failed to match name %s in %s to %s", sqioId(inIO), argv[1], argv[0]) ;
      if (id != sqioId(inIO)) free (id) ;
      if (inIO->seqLen != gdb->seqLen[seq])
	die ("length mismatch for seq %s: %d in %s, %d in %s",
	     sqioId(inIO), gdb->seqLen[seq], argv[0], inIO->seqLen, argv[1]) ;
      U64   nNew = 0, nOld = 0 ;
      char *oldSeq = sqioSeq(inIO) ;
      // printf ("sequence %d %s length %d\n", seq, sqioId(inIO), (int) inIO->seqLen) ;
      for (i = seqStart[seq] ; i < arrayMax(at) && arrp(at,i,TanLine)->seq == seq ; ++i)
	{ TanLine *t = arrp(at,i,TanLine) ;
	  // printf ("  copying %d = %d - %d + %d from %d to %d\n", (int)(t->start - nOld + t->unit),
	  // (int)t->start, (int)nOld, (int)t->unit, (int)nOld, (int)nNew) ;
	  if (t->start < nOld)
	    { t->start = nOld ;
	      if (t->start + t->unit > t->end)
		warn ("buried repeat %s size %d ending %lld when last ended at %lld",
		      sqioId(inIO), t->unit, (long long) t->end, (long long) nOld) ;
	      continue ; // just ignore this 
	    }
	  memcpy (tacoSeq+nNew, oldSeq+nOld, t->start - nOld + t->unit) ;
	  nNew += t->start - nOld + t->unit ;
	  nOld = t->end ;
	}
      memcpy (tacoSeq+nNew, oldSeq+nOld, gdb->seqLen[seq] - nOld) ;
      nNew += gdb->seqLen[seq] - nOld ;
      seqIOwrite (outIO, sqioId(inIO), inIO->descLen ? sqioDesc(inIO) : 0, nNew, tacoSeq, 0) ;
    }
  newFree (tacoSeq, maxSeqLen, char) ;

  newFree (seqStart, gdb->nSeq, I32) ;
  arrayDestroy (at) ;
  gdbDestroy (gdb) ;
  seqIOclose (inIO) ;
  seqIOclose (outIO) ;
  if (outFileName != argv[1]) newFree (outFileName, strlen(argv[1]) + 10, char) ;
  
  return 0 ;
}
