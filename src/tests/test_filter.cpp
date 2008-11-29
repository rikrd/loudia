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

#include "filter.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <iostream>

using namespace std;

int main() {
  int nsamples = 1024;
  int nchannels = 10;

  // Initialize out_samples space
  int out_nchannels;
  int out_nsamples;
  Real* out_samples = new Real[nsamples*nchannels];

  // Initialize constant input
  Real* in_samples = new Real[nsamples*nchannels];  
  for (int j = 0; j < nsamples; j++) {
    for (int k = 0; k < nchannels; k++) {
      in_samples[j*nchannels+k] = 3.0;
    }
  }
  
  // Initialize coefficients
  int n_acoeffs = 1;
  int n_bcoeffs = 10;
  Real* acoeffs = new Real[n_acoeffs*nchannels];
  Real* bcoeffs = new Real[n_bcoeffs*nchannels];
  
  for (int j = 0; j < n_acoeffs; j++) {
    for (int k = 0; k < nchannels; k++) {
      acoeffs[j*nchannels+k] = 1.0;
    }
  }

  for (int j = 0; j < n_bcoeffs; j++) {
    for (int k = 0; k < nchannels; k++) {
      bcoeffs[j*nchannels+k] = 1/Real(n_bcoeffs);
    }
  }
  


  Filter flt(acoeffs, n_acoeffs, bcoeffs, n_bcoeffs, 44100, nchannels);
  flt.setup();

  for (int i=0; i<100; i++) {   
    flt.process(in_samples, 1024, 30, &out_samples, &out_nsamples, &out_nchannels);
  }

  delete[] in_samples;

  return 0;
}
