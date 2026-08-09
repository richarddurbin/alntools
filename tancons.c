
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

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  static char *usage = "Usage: tancons [-o <output_file>] [-u <unit_size>] [-s <sequence_file>] [-c <seq>:<start>-<end>] <file[.1ano]>" ;
  if (!argc) { fprintf (stderr, "%s\n", usage) ; exit(1) ; }

  FILE *outFile = stdout ;
  int   unit    = 0 ;
  char *seqFile = 0 ;
  char *seqId   = 0 ;
  I64   startCoord = 0, endCoord = 0 ;

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
        { seqId = *++argv ; parseCoord (seqId, &startCoord, &endCoord) ; --argc ; }
      else die ("unknown option %s\n%s", *argv, usage) ;
      --argc ; ++argv ;
    }
  if (!seqId && unit <= 0) 
    die ("must specify either unit size with -u or sequence coordinates with -c\n%s", usage) ;

  if (argc != 1) die (usage) ;
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (*argv, schema, "ano", 1) ;
  if (!of) die ("failed to open .1ano file %s", *argv) ;
  oneSchemaDestroy (schema) ; // no longer needed
  Gdb *gdb = readGdb (of, 1, stderr) ;
    
  SeqIO *sio ;
  if (seqFile) sio = seqIOopenRead (seqFile, dna2indexConv, "r") ;
  else sio = seqIOopenRead (gdb->seqFileName, dna2indexConv, "r") ;
  if (!sio) die ("failed to open sequence file %s", seqFile) ;
  
  U32 seqNum = 0 ;
  if (seqId)
    { if (!dictFind (gdb->seqDict, seqId, &seqNum))
        die ("failed to find %s in sequence names for %s", seqId, *argv) ;
      if (endCoord && endCoord > gdb->seqLen[seqNum])
        die ("end coord %lld exceeds length of sequence %s (%lld)", endCoord, seqId, gdb->seqLen[seqNum]) ;
    }
  else // find the longest tandem repeat of the given unit size
    { I64 maxLength = 0 ;
      if (!oneGoto (of,'M', 1) || !oneReadLine (of)) 
        die ("failed to find M line in .1ano file") ;
      while (of->lineType == 'M')
        { int seq = oneInt(of,0) ;
          I64 start = oneInt(of,1) ;
          I64 len = oneInt(of,2) - oneInt(of,1) ;
          int u = 0 ;
          while (oneReadLine(of) && of->lineType != 'M')
            if (of->lineType == 'L') u = atoi(oneString(of)) ;
          if (u == unit && len > maxLength) 
            { maxLength = len ;
              seqNum = seq ;
              seqId = dictName(gdb->seqDict, seqNum) ;
              startCoord = start ;
            }
        }
      fprintf (stderr,"longest tandem repeat of unit %d is %lld bases in sequence %s starting at %lld\n", 
               unit, (long long)maxLength, seqId, (long long)startCoord) ;
    }
  
  I64 nParse, nStart = 0, *starts ;
  if (!oneGoto (of, 'M', 1)) die ("can't locate to first M line in .1ano file") ;
  while (!nStart && oneReadLine (of))
    if (of->lineType == 'M' && oneInt(of,0) == seqNum && oneInt(of,1) == startCoord)
      while (oneReadLine (of))
        { if (of->lineType == 'L')
            { if (!unit) unit = atoi (oneString(of)) ;
              else if (unit != atoi (oneString(of)))
                die ("unit mismatch %s != %d\n", oneString(of), unit) ;
            }
          else if (of->lineType == 'P' && unit)
            { nParse = oneLen(of) ;
              I64 *parse = oneIntList(of) ;
              starts = new (nParse, I64) ;
              for (int i = 0 ; i < nParse-1 ; ++i)
                if (parse[i+1] - parse[i] == unit) starts[nStart++] = startCoord + parse[i] ;
              break ; // we are done with the longest match, no need to read further
            }
          else if (of->lineType == 'M') die ("missing L or P line in .1ano file") ;
        }
  fprintf (stderr,"found %lld starts of exact unit size, ", (long long)nStart) ;
 
  // now find the starts in the sequence file and build the consensus
  int i, j ;
  int counts[unit][5], sumCount = 0 ; // A,C,G,T,N
  memset (counts, 0, unit*5*sizeof(int)) ;
  char seq[unit] ;
  while (seqIOread(sio))
    if (!strcmp(sqioId(sio), seqId))
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
  fprintf (stderr,"with average_match %.1f%% in sequence file %s\n", 
           100.0*sumCount/(nStart*unit), seqFile ? seqFile : gdb->seqFileName) ; 
  fprintf (outFile, ">cons%d\n%s\n", unit, seq) ;
  
  gdbDestroy (gdb) ;
  oneFileClose (of) ;
  newFree (starts, nParse, I64) ;
}