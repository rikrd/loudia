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

BandPass::BandPass( int order, Real startFrequency, Real stopFrequency, FilterType filterType, Real passRipple, Real stopAttenuation ) : 
  _order( order ),
  _startFrequency( startFrequency ),
  _stopFrequency( stopFrequency ),
  _passRipple( passRipple ),
  _stopAttenuation( stopAttenuation ),
  _channels( 1 ),
  _filter( _channels ),
  _filterType( filterType ) 
{
  DEBUG("BANDPASS: Constructor order: " << _order 
        << ", startFrequency: " << _startFrequency
        << ", stopFrequency: " << _stopFrequency
        << ", passRipple: " << _passRipple
        << ", stopAttenuation: " << _stopAttenuation );

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
    chebyshev1(_order, _passRipple, _channels, &zeros, &poles, &gain);
    break;

  case CHEBYSHEVII:
    chebyshev2(_order, _stopAttenuation, _channels, &zeros, &poles, &gain);
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
  Real warped = 2.0 * fs * tan( M_PI * _startFrequency / fs );
  Real warpedStop = 2.0 * fs * tan( M_PI * _stopFrequency / fs );

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

void BandPass::a(MatrixXR* a) const{
  _filter.a(a);
}

void BandPass::b(MatrixXR* b) const{
  _filter.b(b);
}

void BandPass::process(const MatrixXR& samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void BandPass::reset(){
  // Initial values
  _filter.reset();
}

int BandPass::order() const{
  return _order;
}

void BandPass::setOrder( int order ){
  _order = order;
  setup();
}

Real BandPass::startFrequency() const{
  return _startFrequency;
}
  
void BandPass::setStartFrequency( Real frequency ){
  _startFrequency = frequency;
  setup();
}

Real BandPass::stopFrequency() const{
  return _stopFrequency;
}
  
void BandPass::setStopFrequency( Real frequency ){
  _stopFrequency = frequency;
  setup();
}

FilterType BandPass::filterType() const{
  return _filterType;
}

void BandPass::setFilterType( FilterType type ){
  _filterType = type;
  setup();
}

Real BandPass::passRipple() const{
  return _passRipple;
}

void BandPass::setPassRipple( Real rippleDB ){
  _passRipple = rippleDB;
  setup();
}

Real BandPass::stopAttenuation() const{
  return _stopAttenuation;
}

void BandPass::setStopAttenuation( Real attenuationDB ){
  _stopAttenuation = attenuationDB;
  setup();
}
