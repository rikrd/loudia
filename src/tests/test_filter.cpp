/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com>
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

#include "Filter.h"

#include "Typedefs.h"
#include "Debug.h"

#include <Eigen/Core>
#include <iostream>

using namespace std;

int main() {
  int nsamples = 1024;
  int nchannels = 10;

  // TODO: test with smaller input samples than coefficients of a filter

  // Initialize out_samples space
  int out_nchannels;
  int out_nsamples;
  MatrixXR out_samples(nsamples, nchannels);

  // Initialize constant input
  MatrixXR in_samples(nsamples, nchannels);  
  for (int j = 0; j < nsamples; j++) {
    for (int k = 0; k < nchannels; k++) {
      in_samples(j, k) = 3.0;
    }
  }
  
  // Initialize coefficients
  int n_acoeffs = 1;
  int n_bcoeffs = 10;
  MatrixXR acoeffs(n_acoeffs, nchannels);
  MatrixXR bcoeffs(n_bcoeffs, nchannels);
  
  for (int j = 0; j < n_acoeffs; j++) {
    for (int k = 0; k < nchannels; k++) {
      acoeffs(j, k) = 1.0;
    }
  }

  for (int j = 0; j < n_bcoeffs; j++) {
    for (int k = 0; k < nchannels; k++) {
      bcoeffs(j, k) = 1/Real(n_bcoeffs);
    }
  }
  
  Filter flt(bcoeffs, acoeffs, nchannels);
  
  for (int i=0; i<100; i++) {   
    flt.process(in_samples, &out_samples);
  }

  return 0;
}
