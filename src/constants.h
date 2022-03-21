#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PI 3.14159265

typedef struct varnames {
  char var[100];
} NAMES;

typedef struct inputvars {
  double *frame;
  double *stats;
  double *vars;
} INVARS;

typedef struct bpix {
  int *xpix;
  int *ypix;
  int blength;
  int xlim1;
  int xlim2;
  int ylim1;
  int ylim2;
} BOUND;

typedef struct apix {
  int *xpix;
  int *ypix;
  int *intensity;
  int npix;
  int width;
  int height;
  int *image;
  int *mask;
  int lev0num;
  int lev1num;
  int lev2num;
  double *lev0Pix;
  double *lev1Pix;
  double *lev2Pix;
  double *cooc01;
  double *cooc12;
  double *cooc02;
} AREA;

#endif
