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

#include "Typedefs.h"
#include "Debug.h"

#include <vector>

#include "BandPass.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

BandPass::BandPass( int order, Real freq, Real freqStop, FilterType filterType, Real ripplePass, Real rippleStop, int channels ) : 
  _order(order),
  _freq(freq),
  _freqStop(freqStop),
  _ripplePass(ripplePass),
  _rippleStop(rippleStop),
  _channels(channels),
  _filter(channels),
  _filterType(filterType) 
{
  DEBUG("BANDPASS: Constructor order: " << _order 
        << ", freq: " << _freq
        << ", freqStop: " << _freqStop
        << ", ripplePass: " << _ripplePass
        << ", rippleStop: " << _rippleStop );

  if ( order < 1 ) {
    // Throw an exception
  }
  
  setup();
  
  DEBUG("BANDPASS: Constructed");
}

void BandPass::setup(){
  DEBUG("BANDPASS: Setting up...");

  DEBUG("BANDPASS: Getting zpk");  
  // Get the lowpass z, p, k
  MatrixXC zeros, poles;
  Real gain;

  switch( _filterType ){
  case CHEBYSHEVI:
    chebyshev1(_order, _ripplePass, _channels, &zeros, &poles, &gain);
    break;

  case CHEBYSHEVII:
    chebyshev2(_order, _rippleStop, _channels, &zeros, &poles, &gain);
    break;

  case BUTTERWORTH:
    butterworth(_order, _channels, &zeros, &poles, &gain);
    break;

  case BESSEL:
    bessel(_order, _channels, &zeros, &poles, &gain);
    break;
  }
  
  DEBUG("BANDPASS: zeros:" << zeros );
  DEBUG("BANDPASS: poles:" << poles );
  DEBUG("BANDPASS: gain:" << gain );
  
  // Convert zpk to ab coeffs
  MatrixXC a;
  MatrixXC b;
  zpkToCoeffs(zeros, poles, gain, &b, &a);

  DEBUG("BANDPASS: Calculated the coeffs");

  // Since we cannot create matrices of Nx0
  // we have created at least one Zero in 0
  if ( zeros == MatrixXC::Zero(zeros.rows(), zeros.cols()) ){
    // Now we must remove the last coefficient from b
    MatrixXC temp = b.block(0, 0, b.rows(), b.cols()-1);
    b = temp;
  }

  // Get the warped critical frequency
  Real fs = 2.0;
  Real warped = 2.0 * fs * tan( M_PI * _freq / fs );
  Real warpedStop = 2.0 * fs * tan( M_PI * _freqStop / fs );

  Real warpedCenter = sqrt(warped * warpedStop);
  Real warpedBandwidth = warpedStop - warped;

  // Warpped coeffs
  MatrixXC wa;
  MatrixXC wb;
  lowPassToBandPass(b, a, warpedCenter, warpedBandwidth, &wb, &wa);

  DEBUG("BANDPASS: Calculated the low pass to band pass");
  
  // Digital coeffs
  MatrixXR da;
  MatrixXR db;
  bilinear(wb, wa, fs, &db, &da);
  
  DEBUG("BANDPASS: setup the coeffs");

  // Set the coefficients to the filter
  _filter.setA( da.transpose() );
  _filter.setB( db.transpose() );
  
  _filter.setup();
  
  DEBUG("BANDPASS: Finished set up...");
}

void BandPass::a(MatrixXR* a) {
  _filter.a(a);
}

void BandPass::b(MatrixXR* b) {
  _filter.b(b);
}

void BandPass::process(MatrixXR samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void BandPass::reset(){
  // Initial values
  _filter.reset();
}
