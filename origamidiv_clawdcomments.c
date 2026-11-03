//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

/* ===========================================================================
 * origamidiv.c  —  [CLAUDE] annotated copy of Veyrat's "subdivide" branch
 *                  (github.com/dveyrat/origami, branch subdivide, 2022)
 *
 * WHAT IT DOES
 *   A version of origamitag.c that computes ORIGAMI morphology (crossing
 *   number: 0 void, 1 wall, 2 filament, 3 halo) for ONE sub-cube of the
 *   simulation instead of the whole box, so a huge simulation (e.g. 2304^3
 *   particles) never has to be held in memory at once.
 *
 *   The box's initial particle grid (np1d per side) is cut into
 *   ndiv x ndiv x ndiv sub-cubes, each npd = np1d/ndiv particles per side.
 *   This program loads one sub-cube PLUS a "buffer" layer nbuf particles
 *   thick on every side (so np2 = npd + 2*nbuf per side), runs the usual
 *   crossing test, then throws the buffer away and saves only the
 *   sub-cube's own particles. The buffer is there so particles near the
 *   sub-cube's edge can still be compared with neighbours that belong to
 *   the next sub-cube over.
 *
 *   Run it once per sub-cube (divid = 0 .. ndiv^3 - 1); the runs are
 *   independent and can go in parallel. origamicomb.c then stitches the
 *   pieces back into one full-box tag file.
 *
 * HOW TO RUN
 *   ./origamidiv parameterfile ndiv divid
 *     ndiv   number of sub-cubes per side (must divide np1d evenly)
 *     divid  which sub-cube to process. Sub-cube position on the grid:
 *              x-index = divid % ndiv
 *              y-index = (divid / ndiv) % ndiv
 *              z-index = divid / (ndiv*ndiv)
 *   Parameter file: same as upstream plus the new "nbuf" line (see
 *   readfiles.c). Note ndiv here comes from the command line, NOT from the
 *   "ndiv" line in the parameter file (that one is for the VOBOZ tessellation).
 *
 * INPUT   posfile (numfiles must be 0; see KNOWN ISSUES). Format read by
 *         divread() in readfiles.c: 8-byte particle count, then all x,
 *         then all y, then all z as 4-byte floats, sorted by particle ID.
 * OUTPUT  outdir/<taglabel>_tag_<divid>.dat : npd^3 bytes, one morphology
 *         value (0-3) per particle, no header. Order: x fastest, then y,
 *         then z, within the sub-cube.
 *
 * KNOWN ISSUES (found while annotating; marked [ISSUE] below)
 *   - Output name "<label>_tag_<divid>.dat" does not match what
 *     origamicomb.c reads ("<label>tag_<x><y><z>.dat"); files must be renamed.
 *   - The Gadget branch (numfiles > 0) loads the whole box with the
 *     non-subdivided reader, which does not fit the sub-cube layout used here.
 *   - argv[2]/argv[3] are read without checking they were given (crash if missing).
 *   - With nsplit > 1, np2 must be divisible by nsplit or some rows are skipped.
 *
 * WHAT VEYRAT CHANGED
 *   This file is a copy of upstream origamitag.c (v2.0). Every line Veyrat
 *   added or changed ends with  [VEYRAT]  (found with a line-by-line diff
 *   against origamitag.c). Untagged code is identical to upstream,
 *   including the whole crossing-detection algorithm. The changes are:
 *     - read ndiv/divid from the command line; compute npd, np2, np3
 *     - grid size np1d -> np2 (sub-cube + buffer) in arrays, par(), goodmod()
 *     - crossing search distance ng4: np1d/4 -> nbuf
 *     - posread() -> divread(); particle count int -> long
 *     - periodic distance wrap "if" -> "while" (same result in practice)
 *     - drop buffer particles into M before writing; per-sub-cube output
 *       file with no header
 *     - two leftover debugging printf lines
 *
 * Comments starting with [CLAUDE] were added; all code is unchanged.
 * See README.md in this folder for a short C primer.
 * ======================================================================== */

#include <stdio.h>     /* [CLAUDE] file reading/writing and printf */
#include <stdlib.h>    /* [CLAUDE] malloc (memory allocation), exit, abs */
#include "params.h"    /* [CLAUDE] names of the global parameters filled by read_parameters() */

/* [CLAUDE] Shorthand macros (text substitutions done before compiling):
 *   DL          loop d over the 3 dimensions (0=x, 1=y, 2=z)
 *   BF          a huge number, used as a starting min/max
 *   max(A,B)    the larger of A and B
 *   goodmod(A,B) wrap index A back into the range 0..B-1 (periodic wrap),
 *               assuming A is at most one period outside that range */
#define DL for (d=0;d<3;d++)
#define BF 1e30
#define max(A,B) (((A)>(B)) ? (A):(B))
#define goodmod(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))

/* [CLAUDE] Returns 1 if h is negative, else 0. Used to step h through
 * 1, -1, 2, -2, 3, -3, ... (see the crossing loops below). */
int isneg(int h) {
    return (int)(h < 0);
}

/* [CLAUDE] Converts 3-D grid coordinates (i,j,k) on a grid with ng points
 * per side into a single array position. x changes fastest, then y, then z. */
int par(int i, int j, int k, int ng) {
    return i + (j + k*ng)*ng;
}

/* [CLAUDE] Declarations of functions defined in readfiles.c, so the compiler
 * knows their argument types. divread() is Veyrat's sub-cube reader. */
void read_parameters(char *filename);
long divread(char *filename, float ***p, float fact, float boxsize, int np1d, int npd, int ndiv, int divid, int nbuf);  /* [VEYRAT] */
int readgadget(char *filename, float ***p);


int main(int argc, char *argv[]) {
    /* [CLAUDE] ---- Variable declarations (C requires listing variables and
     * their types up front). int = whole number, long = big whole number,
     * float = decimal number, unsigned char = one byte (0-255),
     * char[80] = text of up to 79 characters, "*" = pointer (address of data,
     * used here for arrays whose size is decided at run time). ---- */
    int exitcode;
    long np;            /* [CLAUDE] particle count reported by the reader (total in the file) */  /* [VEYRAT] */
    float **r;          /* [CLAUDE] particle positions: r[i][0..2] = x,y,z of particle i */

    FILE *pos, *tag;
    char *paramfile,tagoutfile[80];
    float xmin,xmax,ymin,ymax,zmin,zmax;

    int s;
    float negb2,b2;     /* [CLAUDE] -boxsize/2 and +boxsize/2, for periodic distance wrapping */
    int ng2,ng4, np2,np3, h, i,i2, j, x,y,z,nhalo,nhalo0,nhalo1,nhalo2,nhaloany, ndiv,divid,npd;  /* [VEYRAT] */
    unsigned char *m,*m0,*m1,*m2,*M, mn,m0n,m1n,m2n; /*Morphology tag */  /* [VEYRAT] */
    /* [CLAUDE] m, m0, m1, m2: one byte per particle (including buffer), one
     *   array per set of 3 orthogonal axes (see the crossing section).
     * M: final morphology for the sub-cube's own particles only (buffer removed). */
    int zstart, ystart, xstart, zend, yend, xend;
    float dx,d1,d2;

    /* [CLAUDE] ---- Read the command line ----
     * argv[] holds the words typed after the program name:
     * argv[1] = parameter file, argv[2] = ndiv, argv[3] = divid. */
    if (argc < 2) {
      printf("Please specify a parameter file\n");
      exit(0);
    } else paramfile = argv[1];
    read_parameters(paramfile);   /* [CLAUDE] fills posfile, outdir, boxsize, np1d, nbuf, ... */

    /* [CLAUDE] sscanf converts text to a number; it returns 1 on success.
     * [ISSUE] argc is only checked for >= 2, so leaving out ndiv or divid crashes here. */
    if (sscanf(argv[2],"%d",&ndiv) != 1) {  /* [VEYRAT] */
        printf("That's no ndiv; try again.\n");  /* [VEYRAT] */
        exit(0);  /* [VEYRAT] */
    }  /* [VEYRAT] */
    if (sscanf(argv[3],"%d",&divid) != 1) {  /* [VEYRAT] */
        printf("That's no divid; try again.\n");  /* [VEYRAT] */
        exit(0);  /* [VEYRAT] */
    }  /* [VEYRAT] */
    /* [CLAUDE] % is the remainder after division; the grid must split evenly. */
    if (np1d%ndiv != 0) {  /* [VEYRAT] */
        printf("ndiv must be factor of np1d.\n");  /* [VEYRAT] */
	exit(0);  /* [VEYRAT] */
    }  /* [VEYRAT] */
    npd = np1d/ndiv;   /* [CLAUDE] particles per side of one sub-cube (no buffer) */  /* [VEYRAT] */

    b2 = boxsize/2.;
    negb2 = -boxsize/2.;

    /* [CLAUDE] Build the output filename, e.g. "/data/n256_tag_5.dat".
     * [ISSUE] origamicomb.c expects "<label>tag_<x><y><z>.dat" instead. */
    sprintf(tagoutfile,"%s%s_tag_%d.dat",outdir,taglabel,divid);  /* [VEYRAT] */

    np2 = npd+nbuf+nbuf;   /* [CLAUDE] particles per side INCLUDING the buffer on both sides */  /* [VEYRAT] */
    np3 = np2*np2*np2;     /* [CLAUDE] total particles loaded (sub-cube + buffer) */  /* [VEYRAT] */

    /* [CLAUDE] ng4 is the maximum lattice distance searched for a crossing.
     * Upstream uses np1d/4 (see Falck+2012 §2.1); here it is nbuf, so the
     * search from a sub-cube particle never runs past the edge of the buffer.
     * ng2 is not used. */
    ng2=nbuf*2;  /* [VEYRAT] */
    ng4=nbuf;  /* [VEYRAT] */

    /* [CLAUDE] ---- Load positions ----
     * divread() fills r with the np3 particles of this sub-cube + buffer,
     * laid out on the np2^3 grid via par(x,y,z,np2).
     * [ISSUE] The readgadget() branch loads the full box in the full-box
     * layout, which does not match the np2 grid used below; use numfiles = 0. */
    if (numfiles > 0) {
      np = readgadget(posfile,&r);
    } else np = divread(posfile,&r,1.,boxsize,np1d,npd,ndiv,divid,nbuf);  /* [VEYRAT] */
    /* [CLAUDE] Note: this prints the total in the file, not the np3 loaded. */
    printf("%ld particles\n",np);fflush(stdout);  /* [VEYRAT] */

    /* [CLAUDE] ---- Allocate and initialise the tag arrays ----
     * malloc(n) reserves n bytes of memory and returns its address.
     * Every tag starts at 1; detecting a crossing multiplies it by a prime
     * (2, 3 or 5), one prime per axis, so later we can tell which axes
     * crossed by checking divisibility. */
    xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
    m = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* [VEYRAT] */
    m0 = (unsigned char *)malloc(np3*sizeof(unsigned char)); /* for the diagonals */  /* [VEYRAT] */
    m1 = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* [VEYRAT] */
    m2 = (unsigned char *)malloc(np3*sizeof(unsigned char));  /* [VEYRAT] */
    M = (unsigned char *)malloc(npd*npd*npd*sizeof(unsigned char));  /* [VEYRAT] */
    for (i=0; i<np3;i++) {  /* [VEYRAT] */
        /* [CLAUDE] Track the position range (only used by a commented-out printf). */
        if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
        if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
        if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];

        m[i] = 1;
        m0[i] = 1;
        m1[i] = 1;
        m2[i] = 1;

	if (i<(npd*npd*npd)) {  /* [VEYRAT] */
	    M[i] = 1;  /* [VEYRAT] */
	}  /* [VEYRAT] */
    }

    /* [CLAUDE] NULL means malloc failed (out of memory). This check comes
     * after m has already been written to, so it would not prevent a crash. */
    if (m==NULL) {
        printf("Morphology array cannot be allocated.\n");
        exit(0);
    }
    //  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
    printf("Calculating ORIGAMI morphology.\n");

    /* [CLAUDE] ==== Crossing detection ====
     * "#pragma omp parallel for" tells the compiler (with -fopenmp) to run the
     * iterations of the next loop on several CPU cores at once. The loaded
     * np2^3 grid is split into nsplit^3 blocks, one loop iteration per block.
     * shared(...) = variables all threads see; private(...) = each thread
     * gets its own copy.
     * [ISSUE] np2/nsplit is whole-number division, so if np2 is not a
     * multiple of nsplit the last few rows are never processed. If those
     * rows are inside the sub-cube (not buffer), their tags stay 0 (void). */
#pragma omp parallel for default(none) shared(np1d, np2, ng4, r, boxsize, negb2, b2, m, m0, m1, m2, nsplit) private (x, y, z, h, dx, i, i2, d1, d2, xstart, xend, ystart, yend, zstart, zend)  /* [VEYRAT] */
    for (s = 0; s < nsplit*nsplit*nsplit; s++) {
        /* [CLAUDE] Turn block number s into block coordinates, then into
         * the grid range [start, end) that this block covers.
         * In C, "/" between whole numbers drops the remainder. */
        zstart = s / (nsplit*nsplit);
        zend = zstart + 1;
        ystart = (s - zstart * nsplit * nsplit)/nsplit;
        yend = ystart + 1;
        xstart = (s - zstart * nsplit * nsplit - ystart * nsplit);
        xend = xstart + 1;
        printf("%d\n",xend);      /* [CLAUDE] leftover debugging output */  /* [VEYRAT] */
        zstart *= (np2/nsplit);  /* [VEYRAT] */
        ystart *= (np2/nsplit);  /* [VEYRAT] */
        xstart *= (np2/nsplit);  /* [VEYRAT] */
        zend *= (np2/nsplit);  /* [VEYRAT] */
        yend *= (np2/nsplit);  /* [VEYRAT] */
        xend *= (np2/nsplit);  /* [VEYRAT] */
	printf("%d\n",xend);          /* [CLAUDE] leftover debugging output */  /* [VEYRAT] */

        /* [CLAUDE] Visit every particle (x,y,z = its position on the
         * INITIAL grid) in this block. Buffer particles are processed too,
         * but their results are discarded at the end. */
        for (x=xstart; x<xend; x++){
            //    printf("%d\n",x);fflush(stdout);
            for (y=ystart; y<yend; y++) {
                for (z=zstart; z<zend; z++) {
                    i = par(x,y,z,np2);   /* [CLAUDE] array index of this particle */  /* [VEYRAT] */

                    /* [CLAUDE] ---- How every test below works ----
                     * Walk along one line of the initial grid, checking the
                     * neighbour h steps away, with h = 1, -1, 2, -2, ...
                     * (h = -h + isneg(h) produces that sequence) until
                     * |h| reaches ng4. A crossing means the neighbour is now
                     * on the WRONG side: it started at +h (ahead) but its
                     * current position is behind (or vice versa). That is
                     * the test "dx*h < 0" (dx and h have opposite signs).
                     * The "while" lines fold dx into -boxsize/2..+boxsize/2
                     * so particles on opposite sides of the periodic box
                     * are measured the short way round.
                     * "break" stops the search at the first crossing found. */

                    /* First just along the Cartesian axes */
                    /* x-direction */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,z,np2);  /* [VEYRAT] */
                        dx = r[i2][0]-r[i][0];
                        while (dx < negb2) dx += boxsize;  /* [VEYRAT] */
                        while (dx > b2) dx -= boxsize;  /* [VEYRAT] */
                        if (dx*h < 0.) {
                            /* [CLAUDE] x crossing: multiply by 2 (only once) */
                            if (m[i] % 2 > 0) {
                                m[i] *= 2;
                            }
                            break;
                        }
                    }
                    /* [CLAUDE] y-direction: same test on the y coordinate; prime 3 */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h, np2),z,np2);  /* [VEYRAT] */
                        dx = r[i2][1]-r[i][1];
                        while (dx < negb2) dx += boxsize;  /* [VEYRAT] */
                        while (dx > b2) dx -= boxsize;  /* [VEYRAT] */
                        if (dx*h < 0.) {
                            /*printf("y:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 3 > 0) {
                                m[i] *= 3;
                            }
                            break;
                        }
                    }
                    /* [CLAUDE] z-direction: prime 5 */
                    for (h=1; abs(h)<ng4; h = -h + isneg(h)) {
                        i2 = par(x,y,goodmod(z+h, np2),np2);  /* [VEYRAT] */
                        dx = r[i2][2]-r[i][2];
                        while (dx < negb2) dx += boxsize;  /* [VEYRAT] */
                        while (dx > b2) dx -= boxsize;  /* [VEYRAT] */
                        if (dx*h < 0.) {
                            /*printf("z:%d %d %d %d %f\n",x,y,z,h,dx);*/
                            if (m[i] % 5 > 0) {
                                m[i] *= 5;
                            }
                            break;
                        }
                    }

                    /* [CLAUDE] ---- 45-degree diagonal axes ----
                     * Falck+2012 found Cartesian axes alone miss about half
                     * the crossings, so three more axis triplets are tested.
                     * Each triplet = one Cartesian axis + the two diagonals
                     * in the plane perpendicular to it. Distance along a
                     * diagonal is the sum (d1+d2) or difference (d1-d2) of
                     * the two coordinate differences.
                     * m0 = diagonals in the y-z plane (pairs with x),
                     * m1 = x-z plane (pairs with y), m2 = x-y plane (pairs with z).
                     * Prime 2 = the "+,+" diagonal, prime 3 = the "+,-" diagonal. */

                    // Now do diagonal directions
                    /* [CLAUDE] x triplet: y-z plane, direction (0,+1,+1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h,np2),goodmod(z+h,np2),np2);  /* [VEYRAT] */
                        d1 = r[i2][1]-r[i][1];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 + d2)*h < 0.) {
                            m0[i] *= 2;
                            break;
                        }
                    }
                    /* [CLAUDE] x triplet: y-z plane, direction (0,+1,-1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(x,goodmod(y+h,np2),goodmod(z-h,np2),np2);  /* [VEYRAT] */
                        d1 = r[i2][1]-r[i][1];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 - d2)*h < 0.) {
                            m0[i] *= 3;
                            break;
                        }
                    }
                    // y
                    /* [CLAUDE] y triplet: x-z plane, direction (+1,0,+1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,goodmod(z+h,np2),np2);  /* [VEYRAT] */
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 + d2)*h < 0.) {
                            m1[i] *= 2;
                            break;
                        }
                    }
                    /* [CLAUDE] y triplet: x-z plane, direction (+1,0,-1) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),y,goodmod(z-h,np2),np2);  /* [VEYRAT] */
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][2]-r[i][2];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 - d2)*h < 0.) {
                            m1[i] *= 3;
                            break;
                        }
                    }
                    // z
                    /* [CLAUDE] z triplet: x-y plane, direction (+1,+1,0) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),goodmod(y+h,np2),z,np2);  /* [VEYRAT] */
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][1]-r[i][1];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 + d2)*h < 0.) {
                            m2[i] *=2;
                            break;
                        }
                    }
                    /* [CLAUDE] z triplet: x-y plane, direction (+1,-1,0) */
                    for (h=1; h<ng4; h = -h + isneg(h)) {
                        i2 = par(goodmod(x+h,np2),goodmod(y-h,np2),z,np2);  /* [VEYRAT] */
                        d1 = r[i2][0]-r[i][0];
                        d2 = r[i2][1]-r[i][1];
                        while (d1 < negb2) d1 += boxsize;  /* [VEYRAT] */
                        while (d1 > b2) d1 -= boxsize;  /* [VEYRAT] */
                        while (d2 < negb2) d2 += boxsize;  /* [VEYRAT] */
                        while (d2 > b2) d2 -= boxsize;  /* [VEYRAT] */
                        if ((d1 - d2)*h < 0.) {
                            m2[i] *= 3;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* [CLAUDE] ==== Count crossed axes and keep only the sub-cube ====
     * For each triplet, count how many of its 3 axes crossed:
     *   (a % p == 0) is 1 if a is divisible by prime p (that axis crossed), else 0.
     *   mn  = Cartesian x,y,z           (m: primes 2,3,5)
     *   m0n = x + the two y-z diagonals (m prime 2, m0 primes 2,3)
     *   m1n = y + the two x-z diagonals (m prime 3, m1 primes 2,3)
     *   m2n = z + the two x-y diagonals (m prime 5, m2 primes 2,3)
     * The morphology is the largest of the four counts (0 void .. 3 halo).
     * The nhalo* counters count halo particles per triplet (just printed).
     * Note: they include buffer particles. */
    nhalo = 0;
    nhalo0 = 0;
    nhalo1 = 0;
    nhalo2 = 0;
    nhaloany = 0;
    j = 0;   /* [CLAUDE] next free slot in the output array M */  /* [VEYRAT] */
    for (i=0;i<np3;i++){  /* [VEYRAT] */
        mn = (m[i]%2 == 0) + (m[i]%3 == 0) + (m[i]%5 == 0);
        m0n = (unsigned char)(m[i]%2 == 0) + (unsigned char)(m0[i]%2 == 0) + (unsigned char)(m0[i]%3 == 0);
        m1n = (unsigned char)(m[i]%3 == 0) + (unsigned char)(m1[i]%2 == 0) + (unsigned char)(m1[i]%3 == 0);
        m2n = (unsigned char)(m[i]%5 == 0) + (unsigned char)(m2[i]%2 == 0) + (unsigned char)(m2[i]%3 == 0);
        m[i] = max(mn,max(m0n,max(m1n,m2n)));
        if (mn == 3) nhalo++;
        if (m0n == 3) nhalo0++;
        if (m1n == 3) nhalo1++;
        if (m2n == 3) nhalo2++;
        if (m[i] == 3) {
            nhaloany++;
        }
	/* [CLAUDE] Recover this particle's grid coordinates (inverse of par())
	 * and shift by nbuf so the sub-cube itself runs 0..npd-1. If any
	 * coordinate is outside that range the particle is in the buffer:
	 * "continue" skips to the next particle without saving it. */
	x = (i%np2)-nbuf;  /* [VEYRAT] */
	if (x < 0 || x >= npd) {  /* [VEYRAT] */
	    continue;  /* [VEYRAT] */
	}  /* [VEYRAT] */
	y = ((i/np2)%np2)-nbuf;  /* [VEYRAT] */
	if (y < 0 || y >= npd) {  /* [VEYRAT] */
	    continue;  /* [VEYRAT] */
	}  /* [VEYRAT] */
	z = (i/(np2*np2))-nbuf;  /* [VEYRAT] */
	if (z < 0 || z >= npd) {  /* [VEYRAT] */
	    continue;  /* [VEYRAT] */
	}  /* [VEYRAT] */
	/* [CLAUDE] Sub-cube particle: save it. Because i runs x-fastest,
	 * M ends up in x-fastest, then y, then z order within the sub-cube. */
	M[j] = m[i];  /* [VEYRAT] */
	++j;  /* [VEYRAT] */
    }
    printf("nhalo=%d,%d,%d,%d,%d\n",nhalo,nhalo0,nhalo1,nhalo2,nhaloany);

    /* Output m */
    /* [CLAUDE] ==== Write the sub-cube's tags ====
     * fopen(..., "w") creates/overwrites the file; fwrite copies npd^3
     * bytes from M into it. Unlike upstream origamitag, no particle-count
     * header is written (that line is commented out); origamicomb.c adds
     * one to the combined file. */
    printf("Writing morphology file to %s.\n",tagoutfile);fflush(stdout);
    tag = fopen(tagoutfile,"w");
    if (tag == NULL) {
        printf("Unable to open %s\n",tagoutfile);
        exit(0);
    }
    //fwrite(&np,1, sizeof(int),tag);  /* [VEYRAT] */
    fwrite(M,npd*npd*npd,sizeof(unsigned char),tag);  /* [VEYRAT] */

    return(1);
}
