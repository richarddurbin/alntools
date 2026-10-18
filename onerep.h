/*  File: onerep.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2026
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 00:14 2026 (rd109)
 * Created: Fri Oct 16 18:10:50 2026 (rd109)
 *-------------------------------------------------------------------
 */

static char *repeatSchemaText = {
  "1 3 def 1 1			schema for repeats\n"
  ".\n"
  "P 3 ano                      genome ANnOtation\n" // build on Gene Myers' .1ano
  "S 3 rep                      REPeat annotation: <1 is genome, <2 is library\n"
  ".\n"
  "O g 0                        a GDB skeleton for a genome\n"
  // NB I am wondering about having one or more field(s) for the species name, and/or
  //    ToLID, and/or assembly accession, and/or taxid.  There is redundancy in these
  //    but that requires external lookup to resolve.  Perhaps all on one line, but it
  //    is conceivable to only have/want information for some of them.
  "G S                          collection of sequences constituting a GDB\n"
  "O S 1 6 STRING               id for a sequence\n"
  "D c 1 6 STRING               chromosome: chrom_name\n"
  "D G 1 3 INT                  gap of given length\n"
  "D C 1 3 INT                  contig of given length\n"
  ".\n"
  "O l 0                        embedded library information\n"
  "G I                          collection of entries each with an identifier and length\n"
  "O I 2 6 STRING 3 INT         identifier, length\n"
  "D Y 1 6 STRING               type of repeat, from some ontology\n"
  "D X 1 3 INT                  terminal repeat length, e.g. LTR of an LTR-RT\n"
  ".\n"
  "O m 1 3 INT                  groups given number of following M objects\n"
  "D i 1 6 STRING               identifier for group, e.g. CEN:1, <gene_name>\n"
  "G M                          \n"
  "O M 3 3 INT 3 INT 3 INT      sequence index, start, end\n" // annotated interval
  "D B 1 3 INT                  bed score (identity as fraction of 1000)\n"
  "D T 2 3 INT 3 INT            tandem repeat: length, library_index (or 0 if none)\n"
  "D P 1 8 INT_LIST             parse of tandem repeat into units\n"
  "D H 2 3 INT 4 REAL           HOR (higher order repeat): period, strength\n"
  "D R 0                        is reverse-complement with respect to library entry\n"
  "D L 3 3 INT 3 INT 3 INT      library match: index, start, end\n"
  "D D 1 3 INT                  TSD length: number of bases to right duplicating to left\n"
  ".\n"
  "P 3 seq                      SEQuence\n"
  "S 3 lib                      repeat LIBrary\n"
  ".\n"
  "O S 1 3 DNA                  sequence of the element\n"
  "D I 1 6 STRING               identifier\n"
  "D Y 1 6 STRING               type of repeat, from some ontology\n"
  "D X 1 3 INT                  terminal repeat length, e.g. LTR of an LTR-RT\n"
  "D U 1 6 STRING               TSD (target site duplication) IUPAC pattern\n"
} ;
