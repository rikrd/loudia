/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

#include "BandFilter.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

BandFilter::BandFilter( int order, Real lowFrequency, Real highFrequency, BandType bandType, FilterType filterType, Real passRipple, Real stopAttenuation ) : 
  _channelCount( 1 ),
  _filter( _channelCount )
{
  DEBUG("BANDFILTER: Constructor order: " << order 
        << ", lowFrequency: " << lowFrequency
        << ", highFrequency: " << highFrequency
        << ", passRipple: " << passRipple
        << ", stopAttenuation: " << stopAttenuation );

  if ( order < 1 ) {
    // Throw an exception
  }

  setOrder( order, false );
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setPassRipple( passRipple, false );
  setStopAttenuation( stopAttenuation, false );
  setFilterType( filterType, false ); 
  setBandType( bandType, false );
  
  setup();
  
  DEBUG("BANDFILTER: Constructed");
}

void BandFilter::setup(){
  DEBUG("BANDFILTER: Setting up...");

  _filter.setChannelCount( _channelCount, false );

  DEBUG("BANDFILTER: Getting zpk");  
  // Get the lowpass z, p, k
  MatrixXC zeros, poles;
  Real gain;

  switch( _filterType ){
  case CHEBYSHEVI:
    chebyshev1(_order, _passRipple, _channelCount, &zeros, &poles, &gain);
    break;

  case CHEBYSHEVII:
    chebyshev2(_order, _stopAttenuation, _channelCount, &zeros, &poles, &gain);
    break;

  case BUTTERWORTH:
    butterworth(_order, _channelCount, &zeros, &poles, &gain);
    break;

  case BESSEL:
    bessel(_order, _channelCount, &zeros, &poles, &gain);
    break;
  }
  
  DEBUG("BANDFILTER: zeros:" << zeros );
  DEBUG("BANDFILTER: poles:" << poles );
  DEBUG("BANDFILTER: gain:" << gain );
  
  // Convert zpk to ab coeffs
  MatrixXC a;
  MatrixXC b;
  zpkToCoeffs(zeros, poles, gain, &b, &a);

  DEBUG("BANDFILTER: Calculated the coeffs");

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

  DEBUG("BANDFILTER: Create the band type filter from the analog prototype");

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

  DEBUG("BANDFILTER: Calculated the low pass to band pass");
  
  // Digital coeffs
  MatrixXR da;
  MatrixXR db;
  bilinear(wb, wa, fs, &db, &da);
  
  DEBUG("BANDFILTER: setup the coeffs");

  // Set the coefficients to the filter
  _filter.setA( da.transpose() );
  _filter.setB( db.transpose() );
  
  _filter.setup();
  
  DEBUG("BANDFILTER: Finished set up...");
}

void BandFilter::a(MatrixXR* a) const{
  _filter.a(a);
}

void BandFilter::b(MatrixXR* b) const{
  _filter.b(b);
}

void BandFilter::process(const MatrixXR& samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void BandFilter::reset(){
  // Initial values
  _filter.reset();
}

int BandFilter::order() const{
  return _order;
}

void BandFilter::setOrder( int order, bool callSetup ){
  _order = order;
  if ( callSetup ) setup();
}

Real BandFilter::lowFrequency() const{
  return _lowFrequency;
}
  
void BandFilter::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real BandFilter::highFrequency() const{
  return _highFrequency;
}
  
void BandFilter::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

BandFilter::FilterType BandFilter::filterType() const{
  return _filterType;
}

void BandFilter::setFilterType( FilterType type, bool callSetup ){
  _filterType = type;
  if ( callSetup ) setup();
}

BandFilter::BandType BandFilter::bandType() const{
  return _bandType;
}

void BandFilter::setBandType( BandType type, bool callSetup ){
  _bandType = type;
  if ( callSetup ) setup();
}

Real BandFilter::passRipple() const{
  return _passRipple;
}

void BandFilter::setPassRipple( Real rippleDB, bool callSetup ){
  _passRipple = rippleDB;
  if ( callSetup ) setup();
}

Real BandFilter::stopAttenuation() const{
  return _stopAttenuation;
}

void BandFilter::setStopAttenuation( Real attenuationDB, bool callSetup ){
  _stopAttenuation = attenuationDB;
  if ( callSetup ) setup();
}
