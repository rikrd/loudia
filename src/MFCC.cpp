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

#include "MFCC.h"
#include "MelBands.h"
#include "DCT.h"

using namespace std;
using namespace Eigen;

MFCC::MFCC(Real lowFrequency, Real highFrequency, int bandCount, Real sampleRate, int fftSize, int coefficientCount, Real minSpectrum, Real power) : 
  _minSpectrum( minSpectrum )
{
  DEBUG("MFCC: Constructor lowFrequency: " << lowFrequency << 
        ", highFrequency: " << highFrequency << 
        ", bandCount: " << bandCount << 
        ", sampleRate: "<< sampleRate << 
        ", fftSize: " << fftSize << 
        ", coefficientCount: " << coefficientCount);

  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setBandCount( bandCount, false );
  setSampleRate( sampleRate, false );
  setFftSize( fftSize, false );
  setCoefficientCount( coefficientCount, false );
  setPower( power, false );
  
  setup();
}

MFCC::~MFCC() {}


void MFCC::setup(){
  // Prepare the buffers
  DEBUG("MFCC: Setting up...");

  _melbands.setFftSize( _fftSize, false );
  _melbands.setSampleRate( _sampleRate, false );
  _melbands.setLowFrequency( _lowFrequency, false );
  _melbands.setHighFrequency( _highFrequency, false );
  _melbands.setBandCount( _bandCount, false );
  _melbands.setup();

  _dct.setInputSize( _bandCount, false );
  _dct.setDctSize( _coefficientCount, false );
  _dct.setup();
  
  reset();
  DEBUG("MFCC: Finished set up...");
}


void MFCC::process(const MatrixXR& spectrum, MatrixXR* mfccCoeffs){
  (*mfccCoeffs).resize(spectrum.rows(), _coefficientCount);
  
  for ( int i = 0; i < spectrum.rows(); i++) {  
    DEBUG("MFCC: Processing Melbands");
    // Process the mel bands on the power of the spectrum
    _melbands.process(spectrum.row(i).cwise().square(), &_bands);
    
    DEBUG("MFCC: Processing Log of bands");
    // Apply a power to the log mel amplitudes as in: http://en.wikipedia.org/wiki/Mel_frequency_cepstral_coefficient
    // V. Tyagi and C. Wellekens
    // On desensitizing the Mel-Cepstrum to spurious spectral components for Robust Speech Recognition
    // in Acoustics, Speech, and Signal Processing, 2005. Proceedings. 
    // IEEE International Conference on, vol. 1, 2005, pp. 529–532.
    _bands = (_bands.cwise() + _minSpectrum).cwise().log() / log(10.0);
    _bands = _bands.cwise().pow(_power);
    
    DEBUG("MFCC: Processing DCT");
    // Process the DCT
    _dct.process(_bands, &_coeffs);

    (*mfccCoeffs).row(i) = _coeffs;
  }

  DEBUG("MFCC: Finished Processing");
}

void MFCC::reset(){
  // Initial values
  _bands = MatrixXR::Zero(1, _bandCount);
  _coeffs = MatrixXR::Zero(1, _coefficientCount);

  _melbands.reset();
  _dct.reset();
}

Real MFCC::lowFrequency() const{
  return _lowFrequency;
}
  
void MFCC::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real MFCC::highFrequency() const{
  return _highFrequency;
}
  
void MFCC::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real MFCC::sampleRate() const{
  return _sampleRate;
}
  
void MFCC::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

int MFCC::coefficientCount() const {
  return _coefficientCount;
}

void MFCC::setCoefficientCount( int count, bool callSetup ) {
  _coefficientCount = count;
  if ( callSetup ) setup();
}

int MFCC::bandCount() const {
  return _bandCount;
}

void MFCC::setBandCount( int count, bool callSetup ) {
  _bandCount = count;
  if ( callSetup ) setup();
}

int MFCC::fftSize() const{
  return _fftSize;
}

void MFCC::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

Real MFCC::power() const{
  return _power;
}
  
void MFCC::setPower( Real factor, bool callSetup ){
  _power = factor;
  if ( callSetup ) setup();
}
