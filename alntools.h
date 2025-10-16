/*  File: alntools.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 16 04:26 2025 (rd109)
 * Created: Thu Oct 16 02:49:56 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "dict.h"
#include "utils.h"
#include "ONElib.h"

static char *alnSchemaText =
  "1 3 def 2 1                 schema for aln and FastGA\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "O s 2 3 INT 6 STRING        length and id for group of sequences = a scaffold\n"
  "G S                         scaffolds (s) group sequence objects (S)\n"
  "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
  "O S 1 3 DNA                 sequence\n"
  "D I 1 6 STRING              identifier of sequence\n"
  ".\n"
  "P 3 aln                     ALIGNMENTS\n"
  "D t 1 3 INT                 trace point spacing in a - global\n"
  ".                           GDB skeleton (may not be presend)\n"
  "O g 0                       groups scaffolds into a GDB skeleton\n"
  "G S                         collection of scaffolds constituting a GDB\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
  "D M 1 8 INT_LIST            mask pair list for a contig\n"
  ".\n"
  "O a 0                       groups A's into a colinear chain\n"
  "G A                         chains (a) group alignment objects (A)\n"
  "D p 2 3 INT 3 INT           spacing in a,b between end of previous alignment and start of next\n"
  ".                           alignment: a_read[beg..end] b_read[beg..end], 0-indexing\n"
  "O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"
  "D L 2 3 INT 3 INT           lengths of sequences a and b\n"
  "D R 0                       flag: reverse-complement sequence b\n"
  "D D 1 3 INT                 differences: number of diffs = substitions + indels\n"
  "D T 1 8 INT_LIST            trace points in b\n"
  "D X 1 8 INT_LIST            number of differences in alignment per trace interval\n"
  "D Q 1 3 INT                 quality: alignment confidence in phred units (currently unused)\n"
  "D E 1 3 INT                 match: number of equal bases (currently unused)\n"
  "D Z 1 6 STRING              cigar string: encodes precise alignment (currently unused)\n"
  "D U 1 3 INT                 putative unit size of a TR alignment (FASTAN)\n"
;

typedef struct {
  int   nSeq, nCtg ;
  DICT *seqDict ;	// names of sequences
  I64  *seqLen ;	// lengths of sequences
  I64  *ctgLen ;	// contig lengths
  int  *child ;         // first contig for each sequence (-1 if sequence has no contigs, all gap)
  int  *parent ;	// parent sequence for each contig
  I64  *offset ;	// offset in parent of each contig
} Gdb ;

