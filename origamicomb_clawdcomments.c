/* ===========================================================================
 * origamicomb.c  —  [CLAUDE] annotated copy of Veyrat's "subdivide" branch
 *                   (github.com/dveyrat/origami, branch subdivide, 2022)
 *
 * WHAT IT DOES
 *   Stitches the per-sub-cube tag files written by origamidiv.c back into a
 *   single full-box tag file, in the same particle order as upstream
 *   origamitag: particle ID = x + y*N + z*N^2 (x fastest), N = num_part.
 *
 *   It does this without loading the whole box: the output is built one
 *   "chunk" at a time. Each chunk is a thin slab of N/df^3 complete z-planes
 *   (df = div_factor), which sits inside one row of sub-cubes along z. For
 *   each chunk it reads the matching piece of every sub-cube in that row,
 *   rearranges the bytes into full-box order, and appends them to the output.
 *
 * HOW TO RUN
 *   ./origamicomb label outdir num_part div_factor
 *     label       the taglabel used for the sub-cube files
 *     outdir      directory holding them (needs a trailing "/")
 *     num_part    particles per side of the full box (np1d)
 *     div_factor  sub-cubes per side (ndiv given to origamidiv)
 *   (No parameter file; this program does not use readfiles.c.)
 *
 * INPUT   outdir/<label>tag_<I><J><K>.dat for every sub-cube, where
 *         I, J, K = the sub-cube's x, y, z position (0 .. div_factor-1).
 *         Each holds (N/df)^3 bytes in x-fastest order (origamidiv output).
 * OUTPUT  outdir/<label>tag.dat : 4-byte particle count, then N^3 bytes of
 *         morphology (0-3) in full-box ID order. This is the format
 *         origamihalo and origamienv read.
 *
 * KNOWN ISSUES (found while annotating; marked [ISSUE] below)
 *   - fwrite(&np3,1,sizeof(int)) is missing its file argument, so this file
 *     does not compile as-is. It should be fwrite(&np3,1,sizeof(int),tag).
 *     There is also no origamicomb target in the Makefile.
 *   - Input names differ from origamidiv's output (<label>_tag_<divid>.dat).
 *     Rename: divid -> I = divid % df, J = (divid/df) % df, K = divid/df^2.
 *     The digits are run together, so this breaks for div_factor >= 10.
 *   - num_part must be divisible by div_factor^3 (e.g. 2304 works for df = 2
 *     or 4, not 3).
 *   - np3 = num_part^3 is a 4-byte int, which overflows above ~1290^3
 *     particles (2304^3 overflows), so the header count is wrong for large
 *     boxes. Upstream origamihalo/env also read this count as an int.
 *   - The output is opened in append mode ("a"): re-running without deleting
 *     the old tag.dat adds a second copy to the end of it.
 *   - No check that the input files actually opened.
 *
 * Comments starting with [CLAUDE] were added; all code is unchanged.
 * See README.md in this folder for a short C primer.
 * ======================================================================== */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

  /* [CLAUDE] sub: a list of open input files, one per sub-cube.
   * FILE ** = "address of a list of file handles". */
  FILE **sub, *tag;
  char *label, *outdir, tagsubfile[80], tagoutfile[80];

  int num_part, div_factor;
  /* [CLAUDE] Size shorthands, with N = num_part, df = div_factor, n = N/df:
   *   np3  = N^3        particles in the full box
   *   df2  = df^2,  df3 = df^3
   *   n1d1 = n          particles per side of one sub-cube
   *   n1d3 = N/df^3     z-planes per output chunk (not used directly)
   *   n2d2 = n^2        particles in one z-plane of one sub-cube
   *   n3d3 = n^3        particles in one sub-cube = particles per output chunk
   *   n3d4 = n^3/df     one sub-cube's share of a chunk times df
   *   n3d5 = n^3/df^2   particles one sub-cube contributes to one chunk
   *                     (= n2d2 * n1d3: its part of the chunk's z-planes) */
  int i,j,k,ijk,l,o,n, np3,df2,df3, n1d1,n1d3,n2d2,n3d3,n3d4,n3d5, irem,idiv;
  int *indices;
  /* [CLAUDE] m1: one chunk's bytes, grouped by sub-cube (as read).
   *          m2: one sub-cube's piece of the chunk, just read from disk. */
  unsigned char *m1, *m2;

  /* [CLAUDE] ---- Command line ----
   * argv[1..4] are the words typed after the program name. sscanf converts
   * text to an integer. No check is made that all four were given. */
  label = argv[1];
  outdir = argv[2];
  sscanf(argv[3],"%d",&num_part);
  sscanf(argv[4],"%d",&div_factor);

  np3 = num_part*num_part*num_part;   /* [CLAUDE] [ISSUE] overflows int for N > ~1290 */
  df2 = div_factor*div_factor;
  df3 = df2*div_factor;

  /* [CLAUDE] Each chunk must be a whole number of z-planes, which needs
   * N divisible by df^3 (% = remainder after division). */
  if (num_part % df3 != 0) {
    printf("Particle grid size must be divisible by number of subdivisions.\n");
    exit(0);
  }

  n1d1 = num_part/div_factor;
  n1d3 = num_part/df3;
  n2d2 = n1d1*n1d1;
  n3d3 = n2d2*n1d1;
  n3d4 = n3d3/div_factor;
  n3d5 = n3d4/div_factor;

  /* [CLAUDE] ==== Build the reordering table ====
   * indices[o] says which byte of m1 belongs at position o of the output
   * chunk. This table is the same for every chunk, so it is built once.
   *
   * Inside m1, the sub-cube with x-position I and y-position J occupies
   * block number (J*df + I), each block n3d5 bytes long, holding that
   * sub-cube's particles in local x-fastest order (lx, ly, lz).
   * In the output the order must be global x fastest, then global y,
   * then z, where global x = I*n + lx and global y = J*n + ly.
   *
   * The nested loops produce exactly that order (innermost = fastest):
   *   n : lx = 0..n-1             (offset n)
   *   l : I  = 0..df-1            (offset I*n3d5, jump to next sub-cube in x)
   *   k : ly = 0..n-1             (offset ly*n)
   *   j : J  = 0..df-1            (offset J*df*n3d5, next row of sub-cubes in y)
   *   i : lz = 0..n1d3-1          (offset lz*n^2, next z-plane)
   * and the sum i+j+k+l+n is where that particle sits in m1.
   * ("while" repeats its block as long as the condition is true;
   *  "x += y" means x = x + y; "++o" means o = o + 1.) */
  indices = (int *)malloc(n3d3*sizeof(int));
  i = 0;
  o = 0;
  while (i<n3d5) {
    j = 0;
    while (j<n3d3) {
      k = 0;
      while (k<n2d2) {
        l = 0;
	    ijk = i+j+k;
	    while (l<n3d4) {
	      for (n=0; n<n1d1; n++) {
            indices[o] = ijk+l+n;
	        ++o;
	      }
	      l += n3d5;
	    }
	    k += n1d1;
      }
      j += n3d4;
    }
    i += n2d2;
  }

  /* [CLAUDE] ==== Open every sub-cube file ====
   * File "<label>tag_<i><j><k>.dat" is stored at list position
   * df^2*i + df*j + k. From how they are read below, i = x-position,
   * j = y-position, k = z-position of the sub-cube.
   * fopen(..., "r") opens for reading. [ISSUE] a missing file is not detected. */
  sub = (FILE **)malloc(df3*sizeof(FILE *));
  for (i=0; i<div_factor; i++) {
    for (j=0; j<div_factor; j++) {
      for (k=0; k<div_factor; k++) {
        sprintf(tagsubfile,"%s%stag_%d%d%d.dat",outdir,label,i,j,k);
	    ijk = (df2*i)+(div_factor*j)+k;
	    sub[ijk] = fopen(tagsubfile,"r");
      }
    }
  }

  m1 = (unsigned char *)malloc(n3d3*sizeof(unsigned char));
  m2 = (unsigned char *)malloc(n3d5*sizeof(unsigned char));
  sprintf(tagoutfile,"%s%stag.dat",outdir,label);
  /* [CLAUDE] "a" = append. [ISSUE] an existing tag.dat is added to, not replaced. */
  tag = fopen(tagoutfile,"a");
  /* [CLAUDE] Header: total particle count, as upstream origamitag writes.
   * [ISSUE] the file argument "tag" is missing here; won't compile. */
  fwrite(&np3,1,sizeof(int));

  /* [CLAUDE] ==== Assemble the output one chunk at a time ====
   * There are df^3 chunks. Chunk i lies in sub-cube z-row idiv = i / df^2
   * (each z-row of sub-cubes contains df^2 chunks). irem is unused. */
  for (i=0; i<df3; i++) {
    irem = i%df2;
    idiv = i/df2;
    o = 0;
    /* [CLAUDE] Read the next n3d5 bytes from each of the df^2 sub-cubes in
     * this z-row: x-position k (inner loop), y-position j (outer loop).
     * Files remember how far they've been read, so each pass picks up the
     * next few z-planes of every sub-cube. The pieces are copied into m1
     * one after another, giving block number j*df + k (as indices[] expects). */
    for (j=0; j<div_factor; j++) {
      for (k=0; k<div_factor; k++) {
	    ijk = (df2*k)+(div_factor*j)+idiv;
	    fread(m2,sizeof(unsigned char),n3d5,sub[ijk]);
	    for (l=0; l<n3d5; l++) {
	      m1[o] = m2[l];
	      ++o;
 	    }
      }
    }
    /* [CLAUDE] Write the chunk in full-box order, one byte at a time,
     * picking bytes from m1 via the reordering table. */
    for (o=0; o<n3d3; o++) {
      fwrite(&m1[indices[o]],1,sizeof(unsigned char),tag);
    }
  }

  /* [CLAUDE] Files are not explicitly closed; the operating system
   * closes them (and flushes the output) when the program ends. */
  return(1);
}
