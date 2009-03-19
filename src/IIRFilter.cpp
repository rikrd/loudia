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

#include "IIRFilter.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

IIRFilter::IIRFilter( int order, Real lowFrequency, Real highFrequency, BandType bandType, FilterType filterType, Real passRipple, Real stopAttenuation ) : 
  _order( order ),
  _lowFrequency( lowFrequency ),
  _highFrequency( highFrequency ),
  _passRipple( passRipple ),
  _stopAttenuation( stopAttenuation ),
  _channels( 1 ),
  _filter( _channels ),
  _filterType( filterType ), 
  _bandType( bandType )
{
  DEBUG("IIRFILTER: Constructor order: " << _order 
        << ", lowFrequency: " << _lowFrequency
        << ", highFrequency: " << _highFrequency
        << ", passRipple: " << _passRipple
        << ", stopAttenuation: " << _stopAttenuation );

  if ( order < 1 ) {
    // Throw an exception
  }
  
  setup();
  
  DEBUG("IIRFILTER: Constructed");
}

void IIRFilter::setup(){
  DEBUG("IIRFILTER: Setting up...");

  DEBUG("IIRFILTER: Getting zpk");  
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
  
  DEBUG("IIRFILTER: zeros:" << zeros );
  DEBUG("IIRFILTER: poles:" << poles );
  DEBUG("IIRFILTER: gain:" << gain );
  
  // Convert zpk to ab coeffs
  MatrixXC a;
  MatrixXC b;
  zpkToCoeffs(zeros, poles, gain, &b, &a);

  DEBUG("IIRFILTER: Calculated the coeffs");

  // Since we cannot create matrices of Nx0
  // we have created at least one Zero in 0
  if ( zeros == MatrixXC::Zero(zeros.rows(), zeros.cols()) ){
    // Now we must remove the last coefficient from b
    MatrixXC temp = b.block(0, 0, b.rows(), b.cols()-1);
    b = temp;
  }

  // Get the warped critical frequency
  Real fs = 2.0;
  Real warped = 2.0 * fs * tan( M_PI * _lowFrequency / fs );
  
  Real warpedStop = 2.0 * fs * tan( M_PI * _highFrequency / fs );
  Real warpedCenter = sqrt(warped * warpedStop);
  Real warpedBandwidth = warpedStop - warped;

  // Warpped coeffs
  MatrixXC wa;
  MatrixXC wb;

  DEBUG("IIRFILTER: Create the band type filter from the analog prototype");

  switch( _bandType ){
  case LOWPASS:
    lowPassToLowPass(b, a, warped, &wb, &wa);
    break;
    
  case HIGHPASS:
    lowPassToHighPass(b, a, warped, &wb, &wa);  
    break;

  case BANDPASS:
    lowPassToBandPass(b, a, warpedCenter, warpedBandwidth, &wb, &wa);
    break;
    
  case BANDSTOP:
    lowPassToBandStop(b, a, warpedCenter, warpedBandwidth, &wb, &wa);
    break;
  }

  DEBUG("IIRFILTER: Calculated the low pass to band pass");
  
  // Digital coeffs
  MatrixXR da;
  MatrixXR db;
  bilinear(wb, wa, fs, &db, &da);
  
  DEBUG("IIRFILTER: setup the coeffs");

  // Set the coefficients to the filter
  _filter.setA( da.transpose() );
  _filter.setB( db.transpose() );
  
  _filter.setup();
  
  DEBUG("IIRFILTER: Finished set up...");
}

void IIRFilter::a(MatrixXR* a) const{
  _filter.a(a);
}

void IIRFilter::b(MatrixXR* b) const{
  _filter.b(b);
}

void IIRFilter::process(const MatrixXR& samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void IIRFilter::reset(){
  // Initial values
  _filter.reset();
}

int IIRFilter::order() const{
  return _order;
}

void IIRFilter::setOrder( int order, bool callSetup ){
  _order = order;
  if ( callSetup ) setup();
}

Real IIRFilter::lowFrequency() const{
  return _lowFrequency;
}
  
void IIRFilter::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real IIRFilter::highFrequency() const{
  return _highFrequency;
}
  
void IIRFilter::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

IIRFilter::FilterType IIRFilter::filterType() const{
  return _filterType;
}

void IIRFilter::setFilterType( FilterType type, bool callSetup ){
  _filterType = type;
  if ( callSetup ) setup();
}

IIRFilter::BandType IIRFilter::bandType() const{
  return _bandType;
}

void IIRFilter::setBandType( BandType type, bool callSetup ){
  _bandType = type;
  if ( callSetup ) setup();
}

Real IIRFilter::passRipple() const{
  return _passRipple;
}

void IIRFilter::setPassRipple( Real rippleDB, bool callSetup ){
  _passRipple = rippleDB;
  if ( callSetup ) setup();
}

Real IIRFilter::stopAttenuation() const{
  return _stopAttenuation;
}

void IIRFilter::setStopAttenuation( Real attenuationDB, bool callSetup ){
  _stopAttenuation = attenuationDB;
  if ( callSetup ) setup();
}
