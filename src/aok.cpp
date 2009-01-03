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

#include "typedefs.h"
#include "debug.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>

#include "aok.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

AOK::AOK(int windowSize, int hopSize, int fftSize, Real normVolume) : _windowSize(windowSize), _hopSize(hopSize), _fftSize(fftSize), _normVolume(normVolume) {

  setup();
}

AOK::~AOK() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void AOK::setup(){
  // Prepare the buffers
  DEBUG("AOK: Setting up...");

  tstep = _hopSize;
  tlag = _windowSize;
  fftlen = _fftSize;
  
  if ( fftlen < (2*tlag) )
    {
      fstep = 2*tlag/fftlen;
      fftlen = 2*tlag;
    }
  else
    {
      fstep = 1;
    }
  
  vol = _normVolume;

  /*	tlag = 64; */		/* total number of rectangular AF lags	*/
  nits = (int) log2((Real) tstep+2);	/* number of gradient steps to take each time	*/
  /*	nits = 2; */
  /*	vol = 2.0; */		/* kernel volume (1.0=Heisenberg limit)	*/
  mu = 0.5;		/* gradient descent factor	*/
  
  forget = 0.001;		/* set no. samples to 0.5 weight on running AF	*/
  nraf = tlag;		/* theta size of rectangular AF	*/
  nrad = tlag;		/* number of radial samples in polar AF	*/
  nphi = tlag;		/* number of angular samples in polar AF */
  outdelay = tlag/2;	/* delay in effective output time in samples	*/
  /* nlag-1 < outdelay < nraf to prevent "echo" effect */

  nlag = tlag + 1;	/* one-sided number of AF lags	*/
  mfft = po2(fftlen);
  slen = 1.42*(nlag-1) + nraf + 3;	/* number of delayed samples to maintain	*/
  
  pi = 3.141592654;
  vol = (2.0*vol*nphi*nrad*nrad)/(pi*tlag);	/* normalize volume	*/
  
  DEBUG("AOK: tlen: " << tlen);
  DEBUG("AOK: nlag: " << nlag);
  DEBUG("AOK: slen: " << slen);
  DEBUG("AOK: xlen: " << xlen);
  DEBUG("AOK: fftlen: " << fftlen);
  DEBUG("AOK: nrad: " << nrad);
  DEBUG("AOK: nraf: " << nraf);
  DEBUG("AOK: nlag: " << nlag);

  kfill((nrad * nphi), 0.0, polafm2);
  kfill((nraf * nlag), 0.0, rectafr);
  kfill((nraf * nlag), 0.0, rectafi);
  kfill(slen, 0.0, xr);
  kfill(slen, 0.0, xi);
  kfill(nphi, 1.0, sigma);
  
  //tlen = xlen + nraf + 2;

  rectamake(nlag, nraf, forget, rar, rai, rarN, raiN);/* make running rect AF parms	*/
  plagmake(nrad, nphi, nlag, plag);
  pthetamake(nrad, nphi, nraf, ptheta, maxrad);	/* make running polar AF parms	*/
  rectrotmake(nraf, nlag, outdelay, rectrotr, rectroti);
  rectopol(nraf, nlag, nrad, nphi, req, pheq);
  
  reset();

  DEBUG("AOK: Finished set up...");
}


void AOK::process(MatrixXC frames, MatrixXR* timeFreqRep){
  outct = 0;

  (*timeFreqRep).resize(frames.rows() / tstep, _fftSize);
  
  // fliplr the frames
  for(int i = 0; i < slen / 2; i++){
    frames.col(i).swap(frames.col(slen - 1 - i));
  }

  for ( int row = 0; row < frames.rows(); row++) {  /*  for each temporal frame of samples  */
    //DEBUG("AOK: Processing, row="<<row);
    // Fill in the input vectors
    //DEBUG("AOK: Processing, setting the xr and xi C arrays, ii="<<ii);
    Eigen::Map<MatrixXR>(xr, 1, frames.cols()) = frames.row(row).real();
    Eigen::Map<MatrixXR>(xi, 1, frames.cols()) = frames.row(row).imag();
    //DEBUG("AOK: Processing, finished setting the xr and xi C arrays");
    
    rectaf(xr, xi, nlag, nraf, rar, rai, rarN, raiN, rectafr, rectafi);

    if ( (row % tstep) == 0 )	/* output t-f slice	*/
      {
        mkmag2((nlag*nraf), rectafr, rectafi, rectafm2);
        polafint(nrad, nphi, nraf, maxrad, nlag, plag, ptheta, rectafm2, polafm2);
        sigupdate(nrad, nphi, nits, vol, mu, maxrad, polafm2, sigma);
        
        for (int i=0; i < nlag-1; i++)	/* for each tau	*/
          {
            tfslicer[i] = 0.0;
            tfslicei[i] = 0.0;
            
            
            for (int j = 0; j < nraf; j++)	/* integrate over theta	*/
              {
                  rtemp = ccmr(rectafr[i*nraf+j], rectafi[i*nraf+j], rectrotr[i*nraf+j], rectroti[i*nraf+j]);
                  rtemp1 = ccmi(rectafr[i*nraf+j], rectafi[i*nraf+j], rectrotr[i*nraf+j], rectroti[i*nraf+j]);
                  
                  rtemp2 = rectkern(i, j, nraf, nphi, req, pheq, sigma);
                  tfslicer[i] = tfslicer[i] + rtemp*rtemp2;
                  tfslicei[i] = tfslicei[i] + rtemp1*rtemp2;
                  /*	fprintf(ofp," %d , %d , %g, %g, %g , %g , %g \n", i,j,rectafr[i*nraf+j],rectafi[i*nraf+j],rtemp,rtemp1,rtemp2); */
              }
          }
        for (int i=nlag-1; i < (fftlen-nlag+2); i++)	/* zero pad for FFT	*/
          {
            tfslicer[i] = 0.0;
            tfslicei[i] = 0.0;
          }
        
        for (int i=(fftlen-nlag+2); i < fftlen; i++)	/* fill in c.c. symmetric half of array	*/
          {
              tfslicer[i] = tfslicer[fftlen - i];
              tfslicei[i] = - tfslicei[fftlen - i];
          }
        
        
        
        fft(fftlen, mfft, tfslicer, tfslicei);
        /*
          DEBUG("AOK: Processing, tfslicer: ");
          for(int b=0; b<fftlen; b++ ){
          DEBUG("AOK: Processing, tfslicer["<<b<<"]="<<tfslicer[b]);
          }
        */
        itemp = fftlen/2 + fstep;
        int col = 0;				/* print output slice	*/
        //DEBUG("AOK: Processing, (*timeFreqRep).shape: " << (*timeFreqRep).rows() << ", " << (*timeFreqRep).cols());
        for (int i=itemp; i < fftlen; i=i+fstep)
          {
            //DEBUG("AOK: Processing, row: " << row << ", col: " << col << ", i: " << i);
            (*timeFreqRep)(outct, col) = tfslicer[i];
            col++;
          }
        for (int i=0; i < itemp; i=i+fstep)
          {
            //DEBUG("AOK: Processing, row: " << row << ", col: " << col << ", i: " << i);
            (*timeFreqRep)(outct, col) = tfslicer[i];
            col++;
          }

        outct = outct + 1;        
      }
  }
}

void AOK::reset(){
  // Initial values
}

int AOK::frameSize() const{
  return 2.42 * _windowSize + 3;
}

int AOK::fftSize() const{
  return _fftSize;
}


/****************************************************************/
/*		fft.c						*/
/*		Douglas L. Jones				*/
/*		University of Illinois at Urbana-Champaign	*/
/*		January 19, 1992				*/
/*								*/
/*   fft: in-place radix-2 DIT DFT of a complex input		*/
/*								*/
/*   input:							*/
/*	n:	length of FFT: must be a power of two		*/
/*	m:	n = 2**m					*/
/*   input/output						*/
/*	x:	Real array of length n with real part of data	*/
/*	y:	Real array of length n with imag part of data	*/
/*								*/
/*   Permission to copy and use this program is granted		*/
/*   as long as this header is included.			*/
/****************************************************************/
void AOK::fft(int n, int m, Real x[], Real y[])
{
	int	i,j,k,n1,n2;
	Real	c,s,e,a,t1,t2;


	j = 0;				/* bit-reverse	*/
	n2 = n/2;
	for (i=1; i < n - 1; i++)
	 {
	  n1 = n2;
	  while ( j >= n1 )
	   {
	    j = j - n1;
	    n1 = n1/2;
	   }
	  j = j + n1;

	  if (i < j)
	   {
	    t1 = x[i];
	    x[i] = x[j];
	    x[j] = t1;
	    t1 = y[i];
	    y[i] = y[j];
	    y[j] = t1;
	   }
	 }


	n1 = 0;				/* FFT	*/
	n2 = 1;

	for (i=0; i < m; i++)
	 {
	  n1 = n2;
	  n2 = n2 + n2;
	  e = -6.283185307179586/n2;
	  a = 0.0;

	  for (j=0; j < n1; j++)
	   {
	    c = cos(a);
	    s = sin(a);
	    a = a + e;

	    for (k=j; k < n; k=k+n2)
	     {
	      t1 = c*x[k+n1] - s*y[k+n1];
	      t2 = s*x[k+n1] + c*y[k+n1];
	      x[k+n1] = x[k] - t1;
	      y[k+n1] = y[k] - t2;
	      x[k] = x[k] + t1;
	      y[k] = y[k] + t2;
	     }
	   }
	 }
	
	    
	return;
}
/*								*/
/*   po2: find the smallest power of two >= input value		*/
/*								*/
int	AOK::po2(int n)
{
	int	m, mm;

	mm = 1;
	m = 0;
	while (mm < n) {
	   ++m;
	   mm = 2*mm;
	}

	return(m);
}
/*								*/
/*   power: compute x^n, x and n positve integers		*/
/*								*/
int AOK::power(int x, int n)
{
	int	i,p;

	p = 1;
	for (i=1; i<=n; ++i)
	     p = p*x;
	return(p);
}
/*								*/
/*   zerofill: set array elements to constant			*/
/*								*/
void AOK::kfill(int len,Real k,Real x[])
{
	int	i;

	for (i=0; i < len; i++)
	  x[i] = k;

	return;
}
/*								*/
/*   cshift: circularly shift an array				*/
/*								*/
void AOK::cshift(int len, Real x[])
{
	int	i;
	Real	rtemp;


	rtemp = x[len-1];

	for (i=len-1; i > 0; i--)
	  x[i] = x[i-1];

	x[0] = rtemp;


	return;
}
/*								*/
/*   cmr: computes real part of x times y			*/
/*								*/
Real	AOK::cmr(Real xr, Real xi, Real yr, Real yi)
{
	Real	rtemp;

	rtemp = xr*yr - xi*yi;

	return(rtemp);
}
/*								*/
/*   cmi: computes imaginary part of x times y			*/
/*								*/
Real	AOK::cmi(Real xr, Real xi, Real yr, Real yi)
{
	Real	rtemp;

	rtemp = xi*yr + xr*yi;

	return(rtemp);
}
/*								*/
/*   ccmr: computes real part of x times y*			*/
/*								*/
Real	AOK::ccmr(Real xr, Real xi, Real yr, Real yi)
{
	Real	rtemp;

	rtemp = xr*yr + xi*yi;

	return(rtemp);
}
/*								*/
/*   ccmi: computes imaginary part of x times y*		*/
/*								*/
Real	AOK::ccmi(Real xr, Real xi, Real yr, Real yi)
{
	Real	rtemp;

	rtemp = xi*yr - xr*yi;

	return(rtemp);
}
/*								*/
/*   rectamake: make vector of poles for rect running AF	*/
/*								*/
void AOK::rectamake(int nlag, int n, Real forget, Real ar[], Real ai[], Real arN[], Real aiN[])
{
	int	i,j;
	Real	trig,decay;
	Real	trigN,decayN;


	trig = 6.283185307/n;
	decay = exp(-forget);

	for (j=0; j < n; j++)
	 {
	  ar[j] = decay*cos(trig*j);
	  ai[j] = decay*sin(trig*j);
	 }

	for (i=0; i < nlag; i++)
	 {
	  trigN = 6.283185307*(n-i);
	  trigN = trigN/n;
	  decayN = exp(-forget*(n-i));

	  for (j=0; j < n; j++)
	   {
	    arN[i*n+j] = decayN*cos(trigN*j);
	    aiN[i*n+j] = decayN*sin(trigN*j);
	   }
	 }


	return;
}
/*								*/
/*   pthetamake: make matrix of theta indices for polar samples	*/
/*								*/
void AOK::pthetamake(int nrad, int nphi, int ntheta, Real ptheta[], int maxrad[])
{
	int	i,j;
	Real	theta,rtemp,deltheta;


	deltheta = 6.283185307/ntheta;

	for (i=0; i < nphi; i++)	/* for all phi ...	*/
	 {
	  maxrad[i] = nrad;

	  for (j = 0; j < nrad; j++)	/* and all radii	*/
	   {
	    theta = - ((4.442882938/nrad)*j)*cos((3.141592654*i)/nphi);
	    if ( theta >= 0.0 )
	       {
	        rtemp = theta/deltheta;
		if ( rtemp > (ntheta/2 - 1) )
		 {
		  rtemp = -1.0;
	          if (j < maxrad[i])  maxrad[i] = j;
		 }
	       }
	      else
	       {
	        rtemp = (theta + 6.283185307)/deltheta;
		if ( rtemp < (ntheta/2 + 1) )
		 {
		  rtemp = -1.0;
	          if (j < maxrad[i])  maxrad[i] = j;
		 }
	       }
		
	    ptheta[i*nrad+j] = rtemp;
	   }
	 }


	return;
}
/*								*/
/*   plagmake: make matrix of lags for polar running AF		*/
/*								*/
void AOK::plagmake(int nrad, int nphi, int nlag, Real plag[])
{
	int	i,j;


	for (i=0; i < nphi; i++)        /* for all phi ...      */
	 {
	  for (j = 0; j < nrad; j++)    /* and all radii        */
	   {
	    plag[i*nrad+j] = mklag(nrad,nphi,nlag,i,j);
	   }
	 }


	return;
}
/*								*/
/*   rectopol: find polar indices corresponding to rect samples	*/
/*								*/
void AOK::rectopol(int nraf, int nlag, int nrad, int nphi, Real req[], Real pheq[])
{
	int	i,j,jt;
	Real	pi,deltau,deltheta,delrad,delphi;


	pi = 3.141592654;

	deltau = sqrt(pi/(nlag-1));
	deltheta = 2.0*sqrt((nlag-1)*pi)/nraf;
	delrad = sqrt(2.0*pi*(nlag-1))/nrad;
	delphi = pi/nphi;

	for (i=0; i < nlag; i++)
	 {
	  for (j=0; j < nraf/2; j++)
	   {
	    req[i*nraf +j] = sqrt(i*i*deltau*deltau + j*j*deltheta*deltheta)/delrad;
	    if ( i == 0 )
	      pheq[i*nraf +j] = 0.0;
	    else pheq[i*nraf +j] = (atan((j*deltheta)/(i*deltau)) + 1.570796327)/delphi;
	   }

	  for (j=0; j < nraf/2; j++)
	   {
	    jt = j - nraf/2;
	    req[i*nraf + nraf/2 + j]  = sqrt(i*i*deltau*deltau + jt*jt*deltheta*deltheta)/delrad;
	    if ( i == 0 )
	      pheq[i*nraf + nraf/2 + j] = 0.0;
	    else pheq[i*nraf + nraf/2 + j] = (atan((jt*deltheta)/(i*deltau)) + 1.570796327)/delphi;
	   }
	 }


	return;
}
/*								*/
/*   rectrotmake: make array of rect AF phase shifts		*/
/*								*/
void AOK::rectrotmake(int nraf, int nlag, Real outdelay, Real rectrotr[], Real rectroti[])
{
	int	i,j;
	Real	twopin;

	twopin = 6.283185307/nraf;


	for (i=0; i < nlag; i++)
	 {
	  for (j=0; j < nraf/2; j++)
	   {
	    rectrotr[i*nraf+j] = cos( (twopin*j)*(outdelay - ((Real) i)/2.0 ) );
	    rectroti[i*nraf+j] = sin( (twopin*j)*(outdelay - ((Real) i)/2.0 ) );
	   }
	  for (j=nraf/2; j < nraf; j++)
	   {
	    rectrotr[i*nraf+j] = cos( (twopin*(j-nraf))*(outdelay - ((Real) i)/2.0 ) );
	    rectroti[i*nraf+j] = sin( (twopin*(j-nraf))*(outdelay - ((Real) i)/2.0 ) );
	   }
	 }


	return;
}
/*								*/
/*   rectaf: generate running AF on rectangular grid;		*/
/*	     negative lags, all DFT frequencies			*/
/*								*/
void AOK::rectaf(Real xr[], Real xi[] , int laglen, int freqlen,  Real alphar[], Real alphai[], Real alpharN[], Real alphaiN[], Real afr[], Real afi[])
{
	int	i,j;
	Real	rtemp,rr,ri,rrN,riN;

	for (i=0; i < laglen; i++)
	 {
	  rr = ccmr(xr[0],xi[0],xr[i],xi[i]);
	  ri = ccmi(xr[0],xi[0],xr[i],xi[i]);

	  rrN = ccmr(xr[freqlen-i],xi[freqlen-i],xr[freqlen],xi[freqlen]);
	  riN = ccmi(xr[freqlen-i],xi[freqlen-i],xr[freqlen],xi[freqlen]);

	  for (j = 0; j < freqlen; j++)
	   {
	    rtemp = cmr(afr[i*freqlen+j],afi[i*freqlen+j],alphar[j],alphai[j]) - cmr(rrN,riN,alpharN[i*freqlen+j],alphaiN[i*freqlen+j]) + rr;
	    afi[i*freqlen + j] = cmi(afr[i*freqlen+j],afi[i*freqlen+j],alphar[j],alphai[j]) - cmi(rrN,riN,alpharN[i*freqlen+j],alphaiN[i*freqlen+j]) + ri;
	    afr[i*freqlen + j] = rtemp;
	   }
	 }


	return;
}
/*								*/
/*   polafint: interpolate AF on polar grid;			*/
/*								*/
void AOK::polafint(int nrad, int nphi, int ntheta, int maxrad[], int nlag, Real plag[], Real ptheta[], Real rectafm2[], Real polafm2[])
{
	int	i,j;
	int	ilag,itheta,itheta1;
	Real	rlag,rtheta,rtemp,rtemp1;


	for (i=0; i < nphi/2; i++)	/* for all phi ...	*/
	 {
	  for (j = 0; j < maxrad[i]; j++)	/* and all radii with |theta| < pi */
	   {
	    ilag = (int) plag[i*nrad+j];
	    rlag = plag[i*nrad+j] - ilag;

	    if ( ilag >= nlag )
	       {
		polafm2[i*nrad+j] = 0.0;
	       }
	      else
	       {
		itheta = (int) ptheta[i*nrad+j];
		rtheta = ptheta[i*nrad+j] - itheta;

		itheta1 = itheta + 1;
		if ( itheta1 >= ntheta )  itheta1 = 0;

		rtemp =  (rectafm2[ilag*ntheta+itheta1] - rectafm2[ilag*ntheta+itheta])*rtheta + rectafm2[ilag*ntheta+itheta];
		rtemp1 =  (rectafm2[(ilag+1)*ntheta+itheta1] - rectafm2[(ilag+1)*ntheta+itheta])*rtheta + rectafm2[(ilag+1)*ntheta+itheta];
		polafm2[i*nrad+j] = (rtemp1-rtemp)*rlag + rtemp;
	       }
	   }
	 }


	for (i=nphi/2; i < nphi; i++)	/* for all phi ...	*/
	 {
	  for (j = 0; j < maxrad[i]; j++)	/* and all radii with |theta| < pi */
	   {
	    ilag = (int) plag[i*nrad+j];
	    rlag = plag[i*nrad+j] - ilag;

	    if ( ilag >= nlag )
	       {
		polafm2[i*nrad+j] = 0.0;
	       }
	      else
	       {
		itheta = (int) ptheta[i*nrad+j];
		rtheta = ptheta[i*nrad+j] - itheta;

		rtemp =  (rectafm2[ilag*ntheta+itheta+1] - rectafm2[ilag*ntheta+itheta])*rtheta + rectafm2[ilag*ntheta+itheta];
		rtemp1 =  (rectafm2[(ilag+1)*ntheta+itheta+1] - rectafm2[(ilag+1)*ntheta+itheta])*rtheta + rectafm2[(ilag+1)*ntheta+itheta];
		polafm2[i*nrad+j] = (rtemp1-rtemp)*rlag + rtemp;
	       }
	   }
	 }


	return;
}
/*								*/
/*   mklag: compute radial sample lag				*/
/*								*/
Real	AOK::mklag(int nrad, int nphi, int nlag, int iphi, int jrad)
{
	Real	delay;

	delay = ((1.414213562*(nlag-1)*jrad)/nrad)*sin((3.141592654*iphi)/nphi);


	return(delay);
}
/*								*/
/*   rectkern: generate kernel samples on rectangular grid	*/
/*								*/
Real	AOK::rectkern(int itau, int itheta, int ntheta, int nphi, Real req[], Real pheq[], Real sigma[])
{
	int	iphi,iphi1;
	Real	kern,tsigma;


	iphi = (int) pheq[itau*ntheta + itheta];
	iphi1 = iphi + 1;
	if (iphi1 > (nphi-1))  iphi1 = 0;
	tsigma = sigma[iphi] + (pheq[itau*ntheta + itheta] - iphi)*(sigma[iphi1]-sigma[iphi]);

/*  Tom Polver says on his machine, exp screws up when the argument of */
/*  the exp function is too large */
	kern = exp(-req[itau*ntheta+itheta]*req[itau*ntheta+itheta]/(tsigma*tsigma));

	return(kern);
}
/*								*/
/*   sigupdate: update RG kernel parameters			*/
/*								*/
void AOK::sigupdate(int nrad, int nphi, int nits, Real vol, Real mu0, int maxrad[], Real polafm2[], Real sigma[])
{
	int	ii,i,j;
	Real	grad[1024],gradsum,gradsum1,tvol,volfac,eec,ee1,ee2,mu;


	for (ii=0; ii < nits; ii++)
	 {
	  gradsum = 0.0;
	  gradsum1 = 0.0;

	  for (i=0; i < nphi; i++)
	   {
	    grad[i] = 0.0;

	    ee1 = exp( - 1.0/(sigma[i]*sigma[i]) );	/* use Kaiser's efficient method */
	    ee2 = 1.0;
	    eec = ee1*ee1;

	    for (j=1; j < maxrad[i]; j++)
	     {
	      ee2 = ee1*ee2;
	      ee1 = eec*ee1;

	      grad[i] = grad[i] + j*j*j*ee2*polafm2[i*nrad+j];
	     }
	    grad[i] = grad[i]/(sigma[i]*sigma[i]*sigma[i]);

	    gradsum = gradsum + grad[i]*grad[i];
	    gradsum1 = gradsum1 + sigma[i]*grad[i];
	   }

	  gradsum1 = 2.0*gradsum1;
	  if ( gradsum < 0.0000001 )  gradsum = 0.0000001;
	  if ( gradsum1 < 0.0000001 )  gradsum1 = 0.0000001;

	  mu = ( sqrt(gradsum1*gradsum1 + 4.0*gradsum*vol*mu0) - gradsum1 ) / ( 2.0*gradsum );


	  tvol = 0.0;

	  for (i=0; i < nphi; i++)
	   {
	    sigma[i] = sigma[i] + mu*grad[i];
	    if (sigma[i] < 0.5)  sigma[i] = 0.5;
/*	    printf("sigma[%d] = %g\n", i,sigma[i]); */
	    tvol = tvol + sigma[i]*sigma[i];
	   }

	  volfac = sqrt(vol/tvol);
	  for (i=0; i < nphi; i++)  sigma[i] = volfac*sigma[i];
	 }


	return;
}
/*								*/
/*   mkmag2: compute squared magnitude of an array		*/
/*								*/
void AOK::mkmag2(int tlen, Real xr[], Real xi[], Real xm2[])
{
	int	i;


	for (i=0; i < tlen; i++)
	 {
	  xm2[i] = xr[i]*xr[i] + xi[i]*xi[i];
	 }

	return;
}
