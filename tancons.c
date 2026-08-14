/*  File: tancons.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: make a consensus for at least one tandem array in a FasTAN .1ano file
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 10 16:08 2026 (rd109)
 * Created: Thu Aug 7 02:48:43 2026 (rd109)
 *-------------------------------------------------------------------
 */

#include "alntools.h"
#include "seqio.h"

static int tanLineCompareCount (const void *a, const void *b)
{ TanLine *ta = (TanLine*)a, *tb = (TanLine*)b ; return tb->count - ta->count ; }

void parseCoord (char *s, I64 *start, I64 *end)
{
  char *text = s ; *start = 0 ; *end = 0 ; 
  while (*s && *s != ':') { if (*s == '\\' && s[1]) ++s ; ++s ; }
  if (!*s) return ; // just the identifier given
  *s++ = 0 ; // terminate the identifier
  if (*s == '-') die ("must have start coord after ':' in %s", text) ;
  *start = strtoll(s, &s, 10) ;
  if (*s++ != '-') die ("must have '-' after start coord in %s", text) ;
  *end = strtoll(s, &s, 10) ;
  if (*s) die ("bad end in %s", text) ;
  if (*end && *start > *end) die ("start %llu > end %llu in %s", *start, *end, text) ;
}

bool readAnoBlock (OneFile *of, Gdb *gdb, TanLine *tl, int minCount) // fills tl
{
  if (of->lineType != 'M') die ("bad call to readAnoLine - at line type %c", of->lineType) ;
  tl->seq = oneInt(of,0) ;
  tl->start = oneInt(of,1) ;
  tl->end = oneInt(of,2) ;
  tl->count = 0 ;
  bool retVal = false ;
  while (oneReadLine(of))
    if (of->lineType == 'L') tl->unit = atoi(oneString(of)) ;
    else if (of->lineType == 'P' && tl->unit >= 20)
      { I64 n = oneLen(of) ;
	I64 *parse = oneIntList(of) ;
	for (int i = 0 ; i < n-1 ; ++i) if (parse[i+1] - parse[i] == tl->unit) ++tl->count ;
	if (tl->count && tl->count >= minCount)
	  { I64 *starts = new (tl->count, I64) ;
	    int j = 0 ;
	    for (int i = 0 ; i < n-1 ; ++i)
	      if (parse[i+1] - parse[i] == tl->unit) starts[j++] = parse[i] ;
	    tl->consensus = (char*)starts ; // abuse consensus field
	    retVal = true ;
	  }
      }
    else if (of->lineType == 'M') break ; // standard end of loop
  return retVal ;
}

void printTanLine (TanLine *tl, Gdb *gdb, bool showStarts)
{
  fprintf (stderr, "seq=%s start=%lld end=%lld unit=%d count=%d\n",
           dictName(gdb->seqDict, tl->seq), (long long)tl->start, (long long)tl->end,
           tl->unit, tl->count) ;
  if (showStarts)
    { fprintf (stderr, "starts:") ;
      for (int i = 0 ; i < tl->count ; ++i) fprintf (stderr, " %lld", (long long)tl->starts[i]) ;
      fprintf (stderr, "\n") ;
    }
}

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: tancons [-o <output_file>] [-u <unit_size>] [-s <sequence_file>] [-c <seq>:<start>-<end>] <file[.1ano]>" ;
  if (!argc) { fprintf (stderr, "%s\n", usage) ; exit(1) ; }

  FILE *outFile  = stdout ;
  int   unit     = 0 ;
  char *seqFile  = 0 ;
  char *seqId    = 0 ;
  int   minCount = 100 ;

  Array atl = arrayCreate (1024, TanLine) ;
  while (argc > 0 && **argv == '-')
    { if (!strcmp (*argv, "-o")) 
        { if (!(outFile = fopen (*++argv, "w"))) die ("failed to open output file %s", *argv) ;
          --argc ;
        }
      else if (!strcmp (*argv, "-u") && argc > 1) 
        { unit = atoi(*++argv) ;
          if (unit <= 0) die ("bad unit size %s", *argv) ;
          --argc ; 
        }
      else if (!strcmp (*argv, "-s") && argc > 1) 
        { seqFile = *++argv ; --argc ; }
      else if (!strcmp (*argv, "-c") && argc > 1) 
        { TanLine *tl = arrayp(atl,0,TanLine) ;
	  seqId = *++argv ; parseCoord (seqId, &tl->start, &tl->end) ;
	  --argc ;
	}
      else die ("unknown option %s\n%s", *argv, usage) ;
      --argc ; ++argv ;
    }
  
  if (argc != 1) die (usage) ;
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "ano", 1) ;
  if (!of) die ("failed to open .1ano file %s", *argv) ;
  oneSchemaDestroy (schema) ; // no longer needed
  Gdb *gdb = readGdb (of, 1, 0) ;
  if (!oneGoto (of,'M', 1) || !oneReadLine (of)) // locate onto first M line
    die ("failed to find M line in .1ano file") ;
    
  SeqIO *sio ;
  if (seqFile) sio = seqIOopenRead (seqFile, dna2indexConv, "r") ;
  else sio = seqIOopenRead (gdb->seqFileName, dna2indexConv, "r") ;
  if (!sio) die ("failed to open sequence file %s", seqFile) ;
  
  if (seqId) // find it and put the full entry into atl[0]
    { TanLine *tl0 = arrp(atl,0,TanLine) ; // this starts by holding start and end
      if (!dictFind (gdb->seqDict, seqId, &tl0->seq))
        die ("failed to find %s in sequence names for %s", seqId, *argv) ;
      TanLine tl ;
      while (of->lineType == 'M')
	if (readAnoBlock(of, gdb, &tl, 0))
	  { if (tl.seq == tl0->seq && tl.start == tl0->start && (!tl0->end || tl.end == tl0->end))
	      { *tl0 = tl ;
		break ; // we are done
	      }
	    else newFree (tl0->consensus, tl0->count, I64) ;
	  }
    }
  else if (unit) // find the unit-sized repeat with most exact-length copies and place in atl[0]
    { int maxCount = 0 ;
      TanLine tl ;
      while (of->lineType == 'M')
	if (readAnoBlock(of, gdb, &tl, maxCount))
	  { if (tl.unit == unit && tl.count > maxCount)
	      { maxCount = tl.count ;
		if (arrayMax(atl))
		  { TanLine *tl0 = arrp(atl,0,TanLine) ;
		    newFree (tl0->consensus, tl0->count,I64) ;
		  }
		array(atl,0,TanLine) = tl ;
	      }
	    else newFree (tl.consensus, tl.count, I64) ;
	  }
      if (!arrayMax(atl)) die ("failed to find any tandem repeat of unit %d in .1ano file", unit) ;
      TanLine *tl0 = arrp(atl,0,TanLine) ;
      fprintf (stderr,"tandem repeat of unit %d with most exact copies (%d) is %s:%lld:%lld\n",
               unit, maxCount, dictName (gdb->seqDict, tl0->seq),
	       (long long)tl0->start, (long long)tl0->end) ;
    }
  else // place in atl all TanLines with at least minCount copies
    { while (of->lineType == 'M')
	if (!readAnoBlock(of, gdb, arrayp(atl,arrayMax(atl),TanLine), minCount)) --arrayMax(atl) ;
      if (!arrayMax(atl))
	die ("failed to find any tandem repeat with at least %d copies in .1ano file", minCount) ;
    }

  if (!arrayMax(atl)) die ("can't find any repeat units with required properties") ;

  // sort atl and build an index from seq -> first entry with that seq
  if (arrayMax(atl) > 1) arraySort (atl, tanLineCompareSeq) ;
  int *atlIndex = new (gdb->nSeq, int) ;
  memset (atlIndex, 0xff, gdb->nSeq*sizeof(int)) ; // sets to -1
  TanLine *tl = arrp(atl,arrayMax(atl),TanLine) ;
  for (int i = arrayMax(atl) ; i-- ;) atlIndex[(--tl)->seq] = i ;
 
  // now find the starts in the sequence file and build the consensus
  int i, j, k ;
  U32 seq ;
  int *accum = new(5*(size_t)I16MAX, int) ;
  while (seqIOread(sio))
    if (dictFind (gdb->seqDict, sqioId(sio), &seq))
      for (i = atlIndex[seq] ; i >= 0 && i < arrayMax(atl) && arrp(atl,i,TanLine)->seq == seq ; ++i)
	{ TanLine *tl = arrp(atl,i,TanLine) ;
	  // printTanLine (tl, gdb, true) ;
	  memset (accum, 0, tl->unit*5*sizeof(int)) ;
	  if (tl->start + tl->starts[tl->count-1] + tl->unit > sio->seqLen)
	    die ("start %d + unit %d exceeds sequence length %d",
		 (int)tl->starts[tl->count-1], tl->unit, (int)sio->seqLen) ;
	  char *s = sqioSeq(sio) ;
	  for (j = 0 ; j < tl->count ; ++j)
	    for (k = 0 ; k < tl->unit ; ++k)
	      ++accum[5*k + s[tl->start + tl->starts[j] + k]] ;
	  newFree (tl->starts, tl->count, I64) ;
	  tl->consensus = new (tl->unit+1, char) ;
	  I32 sum = 0 ;
	  for (k = 0 ; k < tl->unit ; ++k)
	    { int max = 0, maxBase = 4;
	      for (int b = 0 ; b < 4 ; ++b)
		if (accum[5*k + b] > max)
		  { max = accum[5*k + b] ;
		    maxBase = b ;
		  }
	      tl->consensus[k] = index2char[maxBase] ;
	      sum += max ;
	    }
	  tl->consensus[tl->unit] = 0 ;
	  tl->score = (1000 * sum) / (tl->unit*tl->count) ;
	}
    else die ("failed to find %s in GDB directory from .1ano file", sqioId(sio)) ;
  seqIOclose (sio) ;

  // now sort on count and output
  arraySort (atl, tanLineCompareCount) ;
  for (i = 0 ; i < arrayMax(atl) ; ++i)
    { TanLine *tl = arrp(atl,i,TanLine) ;
      fprintf (outFile, ">cons%d#%s:%lld-%lld#%d#%.1f\n%s\n", tl->unit,
	       dictName(gdb->seqDict,tl->seq), tl->start, tl->end, tl->count, tl->score*0.1,
	       tl->consensus) ;
    }
  if (outFile != stdout) fclose (outFile) ;

  // finally free up memory used
  newFree (accum, 5*(size_t)I16MAX, int) ;
  for (i = 0 ; i < arrayMax(atl) ; ++i)
    { TanLine *tl = arrp(atl,i,TanLine) ;
      newFree (tl->consensus, tl->unit+1, char) ;
    }
  arrayDestroy (atl) ;
  newFree (atlIndex, gdb->nSeq, int) ;
  gdbDestroy (gdb) ;
  oneFileClose (of) ;
}

/*************** end of file *****************/
