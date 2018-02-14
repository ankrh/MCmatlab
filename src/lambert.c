/* specfunc/lambert.c
 * 
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* Started with code donated by K. Briggs; added
 * error estimates, GSL foo, and minor tweaks.
 * Some Lambert-ology from
 *  [Corless, Gonnet, Hare, and Jeffrey, "On Lambert's W Function".]
 */

// Shamelessly modified for simpler integration into mcxyz by Anders Hansen, DTU Fotonik, February 2018

/* Halley iteration (eqn. 5.12, Corless et al) */
static double
halley_iteration(
  double x,
  double w_initial,
  unsigned int max_iters
  )
{
  double w = w_initial;
  unsigned int i;

  for(i=0; i<max_iters; i++) {
    double tol;
    const double e = exp(w);
    const double p = w + 1.0;
    double t = w*e - x;
    /* printf("FOO: %20.16g  %20.16g\n", w, t); */

    if (w > 0) {
      t = (t/p)/e;  /* Newton iteration */
    } else {
      t /= e*p - 0.5*(p + 1.0)*t/p;  /* Halley iteration */
    };

    w -= t;

    tol = 10 * DBL_EPSILON * fmax(fabs(w), 1.0/(fabs(p)*e));

    if(fabs(t) < tol) return w;
  }

  /* should never get here */
  return 0;
}


/* series which appears for q near zero;
 * only the argument is different for the different branches
 */
static double
series_eval(double r)
{
  static const double c[12] = {
    -1.0,
     2.331643981597124203363536062168,
    -1.812187885639363490240191647568,
     1.936631114492359755363277457668,
    -2.353551201881614516821543561516,
     3.066858901050631912893148922704,
    -4.175335600258177138854984177460,
     5.858023729874774148815053846119,
    -8.401032217523977370984161688514,
     12.250753501314460424,
    -18.100697012472442755,
     27.029044799010561650
  };
  const double t_8 = c[8] + r*(c[9] + r*(c[10] + r*c[11]));
  const double t_5 = c[5] + r*(c[6] + r*(c[7]  + r*t_8));
  const double t_1 = c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
  return c[0] + r*t_1;
}


double
gsl_sf_lambert_W0(double x)
{
  const double q = x + exp(-1);

  if(x == 0.0) {
    return 0.0;
  }
  else if(q < 0.0) {
    /* Strictly speaking this is an error. But because of the
     * arithmetic operation connecting x and q, I am a little
     * lenient in case of some epsilon overshoot. The following
     * answer is quite accurate in that case. Anyway, we have
     * to return GSL_EDOM.
     */
    return -1.0;
  }
  else if(q == 0.0) {
    return -1.0;
  }
  else if(q < 1.0e-03) {
    /* series near -1/E in sqrt(q) */
    const double r = sqrt(q);
    return series_eval(r);
  }
  else {
    static const unsigned int MAX_ITERS = 10;
    double w;

    if (x < 1.0) {
      /* obtain initial approximation from series near x=0;
       * no need for extra care, since the Halley iteration
       * converges nicely on this branch
       */
      const double p = sqrt(2.0 * exp(1) * q);
      w = -1.0 + p*(1.0 + p*(-1.0/3.0 + p*11.0/72.0)); 
    }
    else {
      /* obtain initial approximation from rough asymptotic */
      w = log(x);
      if(x > 3.0) w -= log(w);
    }

    return halley_iteration(x, w, MAX_ITERS);
  }
}


double
gsl_sf_lambert_Wm1(double x)
{
  if(x > 0.0) {
    return gsl_sf_lambert_W0(x);
  }
  else if(x == 0.0) {
    return 0.0;
  }
  else {
    static const unsigned int MAX_ITERS = 32;
    const double q = x + exp(-1);
    double w;

    if (q < 0.0) {
      /* As in the W0 branch above, return some reasonable answer anyway. */
      return -1.0;
    }

    if(x < -1.0e-6) {
      /* Obtain initial approximation from series about q = 0,
       * as long as we're not very close to x = 0.
       * Use full series and try to bail out if q is too small,
       * since the Halley iteration has bad convergence properties
       * in finite arithmetic for q very small, because the
       * increment alternates and p is near zero.
       */
      const double r = -sqrt(q);
      w = series_eval(r);
      if(q < 3.0e-3) {
        /* this approximation is good enough */
        return w;
      }
    }
    else {
      /* Obtain initial approximation from asymptotic near zero. */
      const double L_1 = log(-x);
      const double L_2 = log(-L_1);
      w = L_1 - L_2 + L_2/L_1;
    }

    return halley_iteration(x, w, MAX_ITERS);
  }
}
