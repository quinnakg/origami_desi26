//-----------------------------------------------------------------------------
// Copyright (c) 2016, Bridget L. Falck, Mark C. Neyrinck, & Nuala McCullagh
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
//-----------------------------------------------------------------------------

/* ===========================================================================
 * readfiles.c  —  [CLAUDE] annotated copy of Veyrat's "subdivide" branch
 *                 (github.com/dveyrat/origami, branch subdivide, 2022)
 *
 * WHAT IT IS
 *   Shared input routines linked into every ORIGAMI/VOBOZ program (it has no
 *   main(); other programs call these functions):
 *     read_parameters()  read the parameter file into global variables
 *     posread()          read a whole-box position file           (upstream)
 *     divread()          read ONE sub-cube + buffer of that file  (Veyrat, NEW)
 *     velread()          read a velocity file                     (upstream)
 *     readgadget() / load_snapshot() / allocate_memory() / reordering()
 *                        read Gadget snapshots and sort by ID     (upstream)
 *
 * WHAT VEYRAT CHANGED (everything else is identical to upstream v2.0)
 *   1. New global parameter nbuf (buffer thickness, in particles) read from
 *      a new section of the parameter file. Because this file is shared,
 *      EVERY program built with it now needs a parameter file that has the
 *      nbuf section, laid out like the branch's parameters_example.txt:
 *          ...
 *          numfiles   0
 *          <blank line>
 *          # Used by origamitag_div
 *          nbuf       64
 *          <blank line>
 *          # Used by vozinit
 *          ...
 *   2. New function divread(), used by origamidiv.c.
 *
 * Comments starting with [CLAUDE] were added; all code is unchanged.
 * See README.md in this folder for a short C primer.
 * ======================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include "params.h"

#define DL for (d=0;d<3;d++) /* Dimension loop */
#define BF 1e30

/* Gadget variables: */
/* [CLAUDE] A "struct" groups several named values into one record.
 * This one mirrors the 256-byte header at the start of each Gadget file. */
struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

int     NumPart, Ngas;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P;   /* [CLAUDE] P: list of particle records, filled by load_snapshot() */

long long *Id;   /* [CLAUDE] Gadget particle IDs (8-byte integers) */

double  Time, Redshift;

/* Global parameters */
/* [CLAUDE] Variables declared outside any function are "global": every
 * program linked with this file can use them (params.h lists them).
 * char * = text; float = decimal number; int = whole number.
 * nbuf is Veyrat's addition. */
int maxline = 80;
char *posfile;
char *outdir;
char *taglabel;
float boxsize;
int np1d;
int nsplit;
int numfiles;
int nbuf;
float buffer;
int numdiv;
float volcut;
int npmin;
char *halolabel;

/* [CLAUDE] ==== read_parameters ====
 * Reads the parameter file in a FIXED order, not by name: each fscanf reads
 * "<word> <value>" and throws the word away (into inbuf). So the lines must
 * appear in exactly the order below, and every one must be present.
 * fgets(...) reads and discards one whole line; groups of three skip the end
 * of the previous line, the blank line, and the "# Used by ..." comment. */
void read_parameters(char *filename) {

  FILE *fd;
  char inbuf[maxline];

  /* [CLAUDE] Reserve 80 characters for each text parameter (paths longer
   * than 79 characters will overflow). */
  posfile = (char *)malloc(maxline*sizeof(char));
  outdir = (char *)malloc(maxline*sizeof(char));
  taglabel = (char *)malloc(maxline*sizeof(char));
  halolabel = (char *)malloc(maxline*sizeof(char));

  fd = fopen(filename,"r");
  if (fd == NULL) {
    printf("Unable to open parameter file %s\n\n",filename);
    exit(0);
  }
  printf("Reading parameter file %s\n",filename);

  fgets(inbuf,maxline,fd);   /* [CLAUDE] skip "# General parameters" */

  fscanf(fd,"%s %s",inbuf,posfile);
  fscanf(fd,"%s %s",inbuf,outdir);
  fscanf(fd,"%s %s",inbuf,taglabel);
  fscanf(fd,"%s %f",inbuf,&boxsize);
  fscanf(fd,"%s %d",inbuf,&np1d);
  fscanf(fd,"%s %d",inbuf,&nsplit);
  fscanf(fd,"%s %d",inbuf,&numfiles);

  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);

  /* [CLAUDE] NEW (Veyrat): buffer thickness in particles for origamidiv.
   * This section sits BEFORE the "# Used by vozinit" section. */
  fscanf(fd,"%s %d",inbuf,&nbuf);

  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);

  fscanf(fd,"%s %f",inbuf,&buffer);
  fscanf(fd,"%s %d",inbuf,&numdiv);

  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);
  fgets(inbuf,maxline,fd);

  fscanf(fd,"%s %f",inbuf,&volcut);
  fscanf(fd,"%s %d",inbuf,&npmin);
  fscanf(fd,"%s %s",inbuf,halolabel);

  fclose(fd);


}


/* Positions */
/* Returns number of particles read */
/* [CLAUDE] ==== posread (upstream) ====
 * Reads the whole box. File layout: 4-byte int count np, then np x values,
 * np y values, np z values (4-byte floats each), sorted by particle ID.
 * Result: p[i][0..2] = x,y,z of particle i, multiplied by fact (unit conversion).
 * "float ***p" lets this function hand the newly created array back to
 * the caller: (*p) is the caller's variable. */
int posread(char *posfile, float ***p, float fact) {

  FILE *pos;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *ptemp;

  //printf("fact=%f\n",fact);
  fflush(stdout);

  pos = fopen(posfile, "r");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,pos); */
  fread(&np,1, sizeof(int),pos); 
  /*fread(&dum,1,4,pos);*/

  /* Allocate the arrays */
  (*p) = (float **)malloc(np*sizeof(float *));
  ptemp = (float *)malloc(np*sizeof(float));

  //printf("np = %d\n",np);

  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) {
    (*p)[i] = (float *)malloc(3*sizeof(float));
    if ((*p)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*p)[i][0] = ptemp[i];
  }
  /*fread(&dum,1,4,pos); 
    fread(&dum,1,4,pos); */
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][1] = ptemp[i];
  /*fread(&dum,1,4,pos);
    fread(&dum,1,4,pos); */
  fread(ptemp,np,4,pos);
  for (i=0; i<np; i++) (*p)[i][2] = ptemp[i];
  /*fread(&dum,1,4,pos); */

  fclose(pos);
  free(ptemp);

  /* Get into physical units */

  for (i=0; i<np; i++) DL (*p)[i][d] *= fact;


  /* Test range -- can comment out */
  /*
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
  */
  return(np);
}

/* Positions in subdivision */
/* Returns number of particles read */
/* [CLAUDE] ==== divread (NEW, Veyrat) ====
 * Reads only the particles of one sub-cube plus a buffer nbuf particles
 * thick around it, selected by their position on the INITIAL grid (worked
 * out from their ID/order in the file, not from where they are now).
 *
 * Arguments:
 *   posfile  position file (layout below)        p      output array
 *   fact     unit-conversion factor              boxsize box side length
 *   np1d     particles per side of full box      npd    per side of a sub-cube
 *   ndiv     sub-cubes per side                  divid  which sub-cube
 *   nbuf     buffer thickness (particles)
 *
 * File layout: 8-byte count np (NOTE: upstream posread reads a 4-byte count,
 * so the same file cannot be used by both), then all x, all y, all z as
 * 4-byte floats, sorted by ID = x + y*np1d + z*np1d^2.
 *
 * Result: p is a grid of np2^3 particles, np2 = npd + 2*nbuf, indexed
 * x-fastest (same as par(x,y,z,np2) in origamidiv.c), with the sub-cube
 * proper at grid positions nbuf .. nbuf+npd-1 on each axis.
 *
 * Periodic wrap: a buffer can stick out past the edge of the box. To handle
 * that, every particle is tried at three copies of the box along each axis
 * (shifted by -1, 0, +1 box lengths, so 27 combinations). Any copy that lands inside
 * the buffered sub-cube is stored, with its position shifted by the same
 * number of box lengths, so neighbouring particles stay next to each other in space.
 *
 * Memory: the file is read in nd3 pieces (np/nd3 values at a time) so the
 * whole coordinate list is never in memory, only the np2^3 kept particles.
 * Returns np, the total in the file (not the number kept). */
long divread(char *posfile, float ***p, float fact, float boxsize, int np1d, int npd, int ndiv, int divid, int nbuf) {

  FILE *pos;
  int dum,d, j,k,l,n, np2,np3,nd3, x,y,z, x2,y2,z2, id;
  long np,i;   /* [CLAUDE] long = 8-byte integer, needed for > 2.1 billion particles */
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *ptemp;

  //printf("fact=%f\n",fact);
  fflush(stdout);

  pos = fopen(posfile, "r");
  if (pos == NULL) {
    printf("Unable to open position file %s\n\n",posfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,pos); */
  /* [CLAUDE] fread(dest, size, count, file) copies count items of size
   * bytes from the file into dest. Here: one 8-byte count. */
  fread(&np,1, sizeof(long),pos);
  /*fread(&dum,1,4,pos);*/

  np2 = npd+nbuf+nbuf;     /* [CLAUDE] particles per side incl. buffer */
  np3 = np2*np2*np2;       /* [CLAUDE] particles kept */
  nd3 = ndiv*ndiv*ndiv;    /* [CLAUDE] number of pieces to read the file in */

  /* Allocate the arrays */
  /* [CLAUDE] p: a list of np3 slots, each will point to that particle's
   * 3 coordinates. ptemp: holding area for one piece of the file. */
  (*p) = (float **)malloc(np3*sizeof(float *));
  ptemp = (float *)malloc((np/nd3)*sizeof(float));

  //printf("np = %d\n",np);

  /* [CLAUDE] ---- Pass 1 of 3: x coordinates ----
   * Loop over the nd3 pieces; read np/nd3 floats into ptemp; then for each
   * particle i in this piece decide where (if anywhere) it goes. */
  for (n=0; n<nd3; n++) {
    fread(ptemp,np/nd3,4,pos);
    for (i=(n*np)/nd3; i<((n+1)*np)/nd3; i++) {
      /* [CLAUDE] Initial-grid coordinates of particle i (from its ID),
       * measured from the corner of this sub-cube, then shifted down by
       * one full box (-np1d) so that adding 0, 1 or 2 boxes below gives
       * the copies at -1, 0, +1 box lengths. */
      x = (i%np1d)-np1d-(npd*(divid%ndiv));
      y = ((i/np1d)%np1d)-np1d-(npd*((divid/ndiv)%ndiv));
      z = (i/(np1d*np1d))-np1d-(npd*(divid/(ndiv*ndiv)));
      /* [CLAUDE] j, k, l pick the z, y, x copy. Adding nbuf converts to
       * coordinates on the buffered grid (0 .. np2-1). Copies falling
       * outside that range are skipped ("continue" = go to next copy). */
      for (j=0; j<3; j++) {
        z2 = z+(np1d*j)+nbuf;
        if (z2<0 || z2>=np2) {
          continue;
        }
        for (k=0; k<3; k++) {
          y2 = y+(np1d*k)+nbuf;
          if (y2<0 || y2>=np2) {
            continue;
          }
          for (l=0; l<3; l++) {
            x2 = x+(np1d*l)+nbuf;
            if (x2<0 || x2>=np2) {
              continue;
            }
            /* [CLAUDE] Slot on the buffered grid (x fastest). In this first
             * pass, reserve space for the particle's 3 coordinates. */
            id = (np2*np2*z2)+(np2*y2)+x2;
            (*p)[id] = (float *)malloc(3*sizeof(float));
            if ((*p)[id] == NULL) {
              printf("Unable to allocate particle array in readfiles!\n");
              fflush(stdout);
              exit(0);
            }
            /* [CLAUDE] Store x, shifted by (l-1) box lengths for the copy
             * used. i-(n*np)/nd3 is i's position within the current piece. */
            (*p)[id][0] = ptemp[i-((n*np)/nd3)]+(l-1)*boxsize;
          }
        }
      }
    }
  }
  /*fread(&dum,1,4,pos); 
    fread(&dum,1,4,pos); */
  /* [CLAUDE] ---- Pass 2 of 3: y coordinates ----
   * Same selection as pass 1 (same particles land in the same slots);
   * the slot already exists, so just store y, shifted by (k-1) boxes. */
  for (n=0; n<nd3; n++){
    fread(ptemp,np/nd3,4,pos);
    for (i=(n*np)/nd3; i<((n+1)*np)/nd3; i++) {
      x = (i%np1d)-np1d-(npd*(divid%ndiv));
      y = ((i/np1d)%np1d)-np1d-(npd*((divid/ndiv)%ndiv));
      z = (i/(np1d*np1d))-np1d-(npd*(divid/(ndiv*ndiv)));
      for (j=0; j<3; j++) {
        z2 = z+(np1d*j)+nbuf;
        if (z2<0 || z2>=np2) {
          continue;
        }
        for (k=0; k<3; k++) {
          y2 = y+(np1d*k)+nbuf;
          if (y2<0 || y2>=np2) {
            continue;
          }
          for (l=0; l<3; l++) {
            x2 = x+(np1d*l)+nbuf;
            if (x2<0 || x2>=np2) {
              continue;
            }
            id = (np2*np2*z2)+(np2*y2)+x2;
            (*p)[id][1] = ptemp[i-((n*np)/nd3)]+(k-1)*boxsize;
          }
        }
      }
    }
  }
  /*fread(&dum,1,4,pos);
    fread(&dum,1,4,pos); */
  /* [CLAUDE] ---- Pass 3 of 3: z coordinates, shifted by (j-1) boxes ---- */
  for (n=0; n<nd3; n++) {
    fread(ptemp,np/nd3,4,pos);
    for (i=(n*np)/nd3; i<((n+1)*np)/nd3; i++) {
      x = (i%np1d)-np1d-(npd*(divid%ndiv));
      y = ((i/np1d)%np1d)-np1d-(npd*((divid/ndiv)%ndiv));
      z = (i/(np1d*np1d))-np1d-(npd*(divid/(ndiv*ndiv)));
      for (j=0; j<3; j++) {
        z2 = z+(np1d*j)+nbuf;
        if (z2<0 || z2>=np2) {
          continue;
        }
        for (k=0; k<3; k++) {
          y2 = y+(np1d*k)+nbuf;
          if (y2<0 || y2>=np2) {
            continue;
          }
          for (l=0; l<3; l++) {
            x2 = x+(np1d*l)+nbuf;
            if (x2<0 || x2>=np2) {
              continue;
            }
            id = (np2*np2*z2)+(np2*y2)+x2;
            (*p)[id][2] = ptemp[i-((n*np)/nd3)]+(j-1)*boxsize;
          }
        }
      }
    }
  }

  /*fread(&dum,1,4,pos); */

  fclose(pos);
  free(ptemp);

  /* Get into physical units */
  /* [CLAUDE] Multiply every coordinate of every kept particle by fact
   * (origamidiv passes 1, i.e. no change). DL = loop over x,y,z. */

  for (i=0; i<np3; i++) DL (*p)[i][d] *= fact;


  /* Test range -- can comment out */
  /*
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np3;i++) {
    if ((*p)[i][0]<xmin) xmin = (*p)[i][0]; if ((*p)[i][0]>xmax) xmax = (*p)[i][0];
    if ((*p)[i][1]<ymin) ymin = (*p)[i][1]; if ((*p)[i][1]>ymax) ymax = (*p)[i][1];
    if ((*p)[i][2]<zmin) zmin = (*p)[i][2]; if ((*p)[i][2]>zmax) zmax = (*p)[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);
  */
  return(np);
}

/* Velocities */
/* Returns number of particles read */
/* [CLAUDE] ==== velread (upstream) ====
 * Same as posread but for velocities (same file layout, 4-byte count).
 * Only used by post-processing, not by the subdivided pipeline. */
int velread(char *velfile, float ***v, float fact) {

  FILE *vel;
  int np,dum,d,i;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float *vtemp;

  vel = fopen(velfile, "r");
  if (vel == NULL) {
    printf("Unable to open velocity file %s\n\n",velfile);
    exit(0);
  }
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Read number of particles */
  /*fread(&dum,1,4,vel);*/ fread(&np,1, sizeof(int),vel); /*fread(&dum,1,4,vel);*/
  printf("np=%d\n",np); fflush(stdout);

  /* Allocate the arrays */
  (*v) = (float **)malloc(np*sizeof(float *));
  vtemp = (float *)malloc(np*sizeof(float));

  printf("about to read velocity array\n"); fflush(stdout);

  /* Fill the arrays */
  /*fread(&dum,1,4,vel);*/
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) {
    (*v)[i] = (float *)malloc(3*sizeof(float));
    if ((*v)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
    (*v)[i][0] = vtemp[i];
  }
  /*fread(&dum,1,4,vel); 
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][1] = vtemp[i];
  /*fread(&dum,1,4,vel);
    fread(&dum,1,4,vel); */
  fread(vtemp,np,4,vel);
  for (i=0; i<np; i++) (*v)[i][2] = vtemp[i];
  /*fread(&dum,1,4,vel); */
  fclose(vel);

  /* Convert from code units into physical units (km/sec) */
  
  //printf("np=%d\n",np); fflush(stdout);
  for (i=0; i<np; i++) DL (*v)[i][d] *= fact;

  /* Test range -- can comment out */
  /*
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if ((*v)[i][0] < xmin) xmin = (*v)[i][0]; if ((*v)[i][0] > xmax) xmax = (*v)[i][0];
    if ((*v)[i][1] < ymin) ymin = (*v)[i][1]; if ((*v)[i][1] > ymax) ymax = (*v)[i][1];
    if ((*v)[i][2] < zmin) zmin = (*v)[i][2]; if ((*v)[i][2] > zmax) zmax = (*v)[i][2];
  }
  printf("vx: %f,%f; vy: %f,%f; vz: %f,%f\n",xmin,xmax, ymin,ymax, zmin,zmax);fflush(stdout);
  */
  return(np);
}


/* [CLAUDE] ==== readgadget (upstream) ====
 * Whole-box reader for Gadget snapshots (used when numfiles > 0): load all
 * files, sort particles by ID, and return pointers to their positions.
 * Not subdivision-aware, so origamidiv.c cannot use it correctly. */
int readgadget(char *filename, float ***pos) {

  int i,d;

  load_snapshot(filename,numfiles);
  printf("NumPart = %d\n",NumPart);

  reordering();  /* call this routine only if your IDs are set properly */

  // Define p from *P:

  /* Allocate the arrays */
  (*pos) = (float **)malloc(NumPart*sizeof(float *));
//  printf("Allocated memory for pos in readgadget\n");

  for (i=0; i<NumPart; i++) {
    (*pos)[i] = (float *)malloc(3*sizeof(float));
    if ((*pos)[i] == NULL) {
      printf("Unable to allocate particle array in readfiles!\n");
      fflush(stdout);
      exit(0);
    }
  }
  printf("Allocated memory for pos in readgadget\n");

  P++; /* start with offset 0 */
  /* [CLAUDE] Point each slot at the position already stored in P (no copy).
   * The 3-float blocks reserved in the loop above are therefore never used. */
  for (i=0; i<NumPart; i++) {
    (*pos)[i] = P[i].Pos;
  }


  return(NumPart);

}



/* (From Volker Springel's Gadget snapshot reader)
 * this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.)
 */
/* [CLAUDE] Each Gadget data block is wrapped in 4-byte markers (a Fortran
 * convention); SKIP reads and discards one marker. For each file
 * "<fname>.<i>" it reads the header, then positions, velocities, IDs,
 * masses and (only for gas) extra gas properties into P and Id. Arrays
 * here start at 1, not 0 (see allocate_memory). */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(long long), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }

  Time= header1.time;
  Redshift= header1.time;
}

/* (From Volker Springel's Gadget snapshot reader)
 * this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(long long))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}

/* (From Volker Springel's Gadget snapshot reader)
 * This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
/* [CLAUDE] Sorts in place by following cycles: put particle i where its
 * ID says it belongs, pick up the particle that was there, put it where
 * ITS ID says, and so on until the cycle closes. Afterwards particle k is
 * at position k, which is what ORIGAMI needs to know its initial grid spot. */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
}

