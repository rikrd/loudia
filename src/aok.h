/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or   
** (at your option) any later version.                                 
**                                                                     
** This program is distributed in the hope that it will be useful,     
** but WITHOUT ANY WARRANTY; without even the implied warranty of      
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
** GNU General Public License for more details.                        
**                                                                     
** You should have received a copy of the GNU General Public License   
** along with this program; if not, write to the Free Software         
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
*/                                                                          

#ifndef AOK_H
#define AOK_H

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class AOK {
protected:
  // Internal parameters
  int _windowSize;
  int _hopSize;
  int _fftSize;
  Real _normVolume;

  // Internal variables
  int itemp;
  int xlen, tlen, nraf,tlag,nlag,mfft,nrad, nphi, nits;
  int slen,fftlen;
  int tstep,fstep,outct;
  
  int maxrad[100000];	/* max radius for each phi s.t. theta < pi */
  
  
  char name[10];
  
  
  Real pi, rtemp, rtemp1, rtemp2, mu, forget;
  Real vol;			/* total kernel volume	*/
  Real outdelay;		/* delay in samples of output time from current sample time	*/
  Real xr[100000],xi[100000];	/* recent data samples	*/
  Real rectafr[70000];		/* real part of rect running AF	*/
  Real rectafi[70000];		/* imag part of rect running AF	*/
  Real rectafm2[70000];	/* rect |AF|^2	*/
  Real polafm2[70000];		/* polar |AF|^2	*/
  Real rectrotr[70000];	/* real part of rect AF phase shift	*/
  Real rectroti[70000];	/* imag part of rect AF phase shift	*/
  Real req[70000];		/* rad corresp. to rect coord	*/
  Real pheq[70000];		/* phi corresp. to rect coord	*/
  Real plag[70000];		/* lag index at polar AF sample points	*/
  Real ptheta[70000];		/* theta index at polar AF sample points	*/
  Real sigma[100000];		/* optimal kernel spreads	*/
  Real rar[100000],rai[100000];	/* poles for running FFTs for rect AF	*/
  Real rarN[70000];		/* poles for running FFTs for rect AF	*/
  Real raiN[70000];
  Real tfslicer[100000];		/* freq slice at current time	*/
  Real tfslicei[100000];

public:
  AOK(int windowSize, int hopSize, int fftSize, Real normVolume = 3.0);

  ~AOK();

  void setup();

  void process(MatrixXC frames, MatrixXR* timeFreqRep);

  void reset();

  int frameSize() const;

  int fftSize() const;

protected:
  void fft(int n, int m, Real x[], Real y[]);

  int po2(int n);

  int power(int x, int n);

  void kfill(int len, Real k, Real x[]);

  void cshift(int len, Real x[]);

  Real cmr(Real xr, Real xi, Real yr, Real yi);

  Real cmi(Real xr, Real xi, Real yr, Real yi);

  Real ccmr(Real xr, Real xi, Real yr, Real yi);

  Real ccmi(Real xr, Real xi, Real yr, Real yi);

  void rectamake(int nlag, int n, Real forget, Real ar[], Real ai[], Real arN[], Real aiN[]);

  void pthetamake(int nrad, int nphi, int ntheta, Real ptheta[], int maxrad[]);

  void plagmake(int nrad, int nphi, int nlag, Real plag[]);

  void rectopol(int nraf, int nlag, int nrad, int nphi, Real req[], Real pheq[]);

  void rectrotmake(int nraf, int nlag, Real outdelay, Real rectrotr[], Real rectroti[]);

  void rectaf(Real xr[], Real xi[] , int laglen, int freqlen,  Real alphar[], Real alphai[], Real alpharN[], Real alphaiN[], Real afr[], Real afi[]);

  void polafint(int nrad, int nphi, int ntheta, int maxrad[], int nlag, Real plag[], Real ptheta[], Real rectafm2[], Real polafm2[]);

  Real mklag(int nrad, int nphi, int nlag, int iphi, int jrad);

  Real rectkern(int itau, int itheta, int ntheta, int nphi, Real req[], Real pheq[], Real sigma[]);

  void sigupdate(int nrad, int nphi, int nits, Real vol, Real mu0, int maxrad[], Real polafm2[], Real sigma[]);

  void mkmag2(int tlen, Real xr[], Real xi[], Real xm2[]);
};

#endif  /* AOK_H */
