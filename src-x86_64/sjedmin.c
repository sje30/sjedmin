#include <R.h>
#include <S.h>			/* for seed_in, seed_out */

/* Use the definitions by Venables & Ripley. */
#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()

/* Sfloat already defined in R.h:55
 * typedef double Sfloat;
 */

void h_amac(Sfloat *t, Sfloat *xd1, Sfloat *p, Sfloat *b, Sfloat *h);

void nnd_n(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b,
	   int *idx, Sfloat *min);

void nnd(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b, int ignore,
	 int *idx, Sfloat *min);

void pairwise_amac(Sfloat *wid, Sfloat *ht, int *numcells,
		   Sfloat *xd1, Sfloat *p, Sfloat *b,
		   Sfloat *xpts, Sfloat *ypts)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */
  
  int i,j, mm, id;
  Sfloat x1, y1, h, dist, u, h1;
  int n = *numcells;

  RANDIN;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = UNIF * (*wid); ypts[i] = UNIF * (*ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    /*Rprintf("select new point for unit %d\n", id);*/
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    do {

      /* generate a new point. */
      xpts[0] = UNIF * (*wid); ypts[0] = UNIF * (*ht);
      u = UNIF;
      h = 1.0;

      for(j=1; j<n; j++) {
	x1 = xpts[j] - xpts[0]; y1 = ypts[j] - ypts[0];
	dist = (sqrt((x1*x1) + (y1*y1)));
	h_amac(&dist, xd1, p, b, &h1);
	/*printf("xxx valu %lf\n", h1);*/
	h *= h1;
	/*Rprintf("t %lf %lf %lf\n", dist, h, h1); */
	/*h = 1;*/
      }
    } while ( h < u);
  }

  RANDOUT;
}


void dminlulbd(Sfloat *wid, Sfloat *ht, int *numcells,
	       Sfloat *pdmin, Sfloat *psd,
	       Sfloat *plower, Sfloat *pupper,
	       Sfloat *xpts, Sfloat *ypts,
	       Sfloat *dmins, int *nrejects)
{
  /* Birth and death version of dmin rule. */
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */
  
  int i,j, mm, id, okay, generate_r;
  Sfloat x1, y1, dist2, u, r, rr, lower, upper;
  int n = *numcells;

  RANDIN;

  lower = *plower; upper = *pupper;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = UNIF * (*wid); ypts[i] = UNIF * (*ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    /*if (i%100==0) Rprintf("it %d:select new point for unit %d\n", i,id);*/
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    do {
      
      /* generate a new point. */
      xpts[0] = UNIF * (*wid); ypts[0] = UNIF * (*ht);
      u = UNIF;
      okay = 1;

      /* genereate a new dmin value. */
      generate_r = 1;
      while (generate_r) {
	r = (*psd * norm_rand()) + *pdmin;
	if ( (r > lower ) && (r < upper)) generate_r = 0; /* okay r value */
      }
      rr = r*r;
      for(j=1; j<n; j++) {
	x1 = xpts[j] - xpts[0]; y1 = ypts[j] - ypts[0];
	dist2 = (x1*x1) + (y1*y1);
	if (dist2 < rr) {
	  okay = 0;
	  break; /* quit for loop and try again. */
	}
      }
    } while ( !okay);

    /* dmins and nrejects not returned here. */
    dmins[id] = r;		/* store r value. */
/*     nrejects[i] = this_cell_rejects; */

  }

  RANDOUT;
}

/* Is this accessible from R? */
void h_amac(Sfloat *t, Sfloat *xd1, Sfloat *p, Sfloat *b, Sfloat *h)
{
  /* Pairwise interaction function for amacrines. */
  if (*t < *xd1)
    *h = 0.0;
  else if (*t <= *p)
    *h = ( pow( ((*t-*xd1)/(*p-*xd1)), *b));
  else
    *h = 1.0;
  /*printf("ret valu %lf\n", *h);*/
}



#define MAXREJECTS 999999

void dminsd(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	    Sfloat *pdmin, Sfloat *psd,
	    Sfloat *xpts, Sfloat *ypts,
	    Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.
   * On exit, we return an array of cells stored in (xs, ys).
   * NJRECTS stores the number of rejected cells.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;


  RANDIN;
  for (i=0; i<*pnumcells; i++) {

    this_dmin = *pdmin + (*psd * norm_rand());
    /*Rprintf("this val %lf\n", this_dmin); */
    dmins[i] = this_dmin;
    looking = 1; this_cell_rejects = 0;
    while (looking) {
      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);
      /*min = 1000; */
      if ( min < this_dmin) {
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminsd error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }

  Rprintf("#rejects %d\tdminsd packing density %.3f\n", num_rejects,
	  ((*pnumcells * PI * *pdmin * *pdmin)/(4*  *pwid * *pht)));
  
  RANDOUT;
}


void dminl(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	   Sfloat *pdmin, Sfloat *psd,
	   Sfloat *xpts, Sfloat *ypts,
	   Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;

  Sfloat min_dmin = 0;		/* Just constrain dmin to be >0 */

  RANDIN;
  /* Just an esitmate since we don't know the actual dmin values */

  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      this_dmin = *pdmin + (*psd * norm_rand());
      if (this_dmin < min_dmin) this_dmin = min_dmin;
      /*Rprintf("this val %lf\n", this_dmin); */

      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);

      if ( min < this_dmin) {
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminl error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }
  
  Rprintf("#rejects %d\tdminl packing density %.3f\n", num_rejects,
	 ((*pnumcells * PI * *pdmin * *pdmin)/(4 * *pwid * *pht)));

  RANDOUT;
}


void dminlul(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	     Sfloat *pdmin, Sfloat *psd,
	     Sfloat *plower, Sfloat *pupper,
	     Sfloat *xpts, Sfloat *ypts,
	     Sfloat *dmins, int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.  Also provide arguments LOWER
   * and UPPER: dmin values outside this range are rejected.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking, generate_dmin;
  int i;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat this_dmin;
  Sfloat lower, upper;



  RANDIN;

  lower = *plower; upper = *pupper;
  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      generate_dmin = 1;
      while (generate_dmin) {
	this_dmin = *pdmin + (*psd * norm_rand());
	if ( (this_dmin > lower) && (this_dmin < upper))
	  generate_dmin = 0;
	/*else printf("dminlul: dmin %f outside range\n", this_dmin);*/
      }

      /*Rprintf("this val %lf\n", this_dmin); */

      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);
      /*min = 1000; */
      if ( (min < this_dmin) ) {
	/* reject cell and try another. */
	num_rejects++;
	this_cell_rejects++;
	if (num_rejects > MAXREJECTS) {
	  /* If we cannot fit more cells in, return and indicate error
	   * by setting first x coordinate to negative value.
	   */
	  printf("dminlul error: num_rejects is too high (%s:%d)\n",
		 __FILE__, __LINE__);
	  xpts[0] = -1;
	  RANDOUT;
	  return;
	}
      }
      else {
	looking = 0;
      }
    }
    /* Accept this cell. */
    dmins[i] = this_dmin;
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }
  
  Rprintf("#rejects %d\tdminl packing density %.3f\n", num_rejects,
	 ((*pnumcells * PI * *pdmin * *pdmin)/( 4 * *pwid * *pht)));

  RANDOUT;
}


void dminacc(Sfloat *pwid, Sfloat *pht, int *pnumcells,
	     Sfloat *acc, Sfloat *dmax, Sfloat *inc,
	     Sfloat *xpts, Sfloat *ypts,
	     int *nrejects)

{
  /* Create NUMCELLS cells distributed in an area of size wid*ht.
   * Minimum distance between cells should be at least dmin.  On exit,
   * we return an array of cells stored in (xs, ys).  NJRECTS stores
   * the number of rejected cells.  This is Lucia's version, where a
   * new dmin value is created every time, regardless of whether the
   * trial cell was accepted or not.
   */
  

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i, dlow, dhigh;
  Sfloat x,y;
  Sfloat min; int idx;
  Sfloat p, plow, phigh;

  RANDIN;
  /* Just an esitmate since we don't know the actual dmin values */

  for (i=0; i<*pnumcells; i++) {

    looking = 1; this_cell_rejects = 0;
    while (looking) {

      x = UNIF * (*pwid); y = UNIF * (*pht);
      
      nnd_n( xpts, ypts, i, x, y, &idx, &min);

      if ( min >  *dmax) {
	looking = 0;
      } else {
	dlow = (int) min; dhigh = dlow+1;
	plow = acc[dlow];
	phigh = acc[dhigh];
	p = plow + ( ((phigh-plow)/ *inc) * (min - dlow));
	if (p > UNIF)
	  looking = 0;		/* accept cell */
	else 
	  this_cell_rejects++;
      }
    }

    /* Accept this cell. */
    xpts[i] = x; ypts[i] = y;
    nrejects[i] = this_cell_rejects;
  }
  
  RANDOUT;
}

void dminacc_bd(Sfloat *wid, Sfloat *ht, int *pnumcells,
		Sfloat *acc, Sfloat *dmax, Sfloat *inc,
		Sfloat *xpts, Sfloat *ypts,
		int *nrejects)

{
  /* Birth and death version of the dminacc function.
   */

  int num_rejects = 0, this_cell_rejects = 0;
  int looking;
  int i, dlow, dhigh;
  int n, mm, id;
  Sfloat x,y, u;
  Sfloat min; int idx;
  Sfloat p, plow, phigh;

  RANDIN;

  n = *pnumcells;

  /* Create initial set of points. */
  for (i=0; i<n; i++) {
    xpts[i] = UNIF * (*wid); ypts[i] = UNIF * (*ht);
  }

  mm = 4 * n;			/* Number of repeats */
  mm *= 10;			/* if starting from Poisson distribution. */
  for (i=0; i<= mm; i++) {
    id = (int)(n*UNIF);	/* pick one point at random. */
    /*if (i%100==0) Rprintf("it %d:select new point for unit %d\n", i,id);*/
    xpts[id] = xpts[0]; ypts[id] = ypts[0];

    looking = 1;
    do {
      
      /* generate a new point. */
      xpts[0] = UNIF * (*wid); ypts[0] = UNIF * (*ht);
      u = UNIF;

      nnd( xpts, ypts, n, xpts[0], ypts[0], 0, &idx, &min);

      if ( min >  *dmax) {
	looking = 0;
      } else {
	dlow = (int) min; dhigh = dlow+1;
	plow = acc[dlow];
	phigh = acc[dhigh];
	p = plow + ( ((phigh-plow)/ *inc) * (min - dlow));
	if (p > UNIF)
	  looking = 0;		/* accept cell */
	else 
	  this_cell_rejects++;
      }
    } while (looking);

    /* Accept this cell. */
    nrejects[id] = this_cell_rejects;
  }
  
  RANDOUT;
}


void nnd_n(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b,
	   int *idx, Sfloat *min)
{
  /* Return the index (IDX) and distance (MIN) of cell closest to point (A,B).
   * The N datapoints are stored in XS[] and YS[].
   */
  if (n==0) {
    *idx = -1; *min = 9999999;
  } else {
    nnd(xs, ys, n, a, b, -1, idx, min);
  }
}

void nnd(Sfloat *xs, Sfloat *ys, int n, Sfloat a, Sfloat b, int ignore,
	int *idx, Sfloat *min)
{
  /* Return the index (IDX) and distance (MIN) of cell closest to point (A,B).
   * The N datapoints are stored in XS[] and YS[].
   * In comparison to nnd_n(), cell number IGNORE is not allowed to be
   * the closest cell.
   */
  int i;
  int index = -1; Sfloat minv = 99999999;

  Sfloat x, y, dx, dy, dist;
  
  for (i=0; i<n; i++) {
    if (i==ignore) {
      continue;
    }
    
    x = xs[i]; y = ys[i];
    dx = x - a; dy = y - b;
    dist = (dx*dx) + (dy*dy);
    if ( dist < minv) {
      index = i;
      minv = dist;
    }
    
  }

  *min = sqrt(minv);
  *idx = index;
}
  
