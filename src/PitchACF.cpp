/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#include "PitchACF.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

PitchACF::PitchACF(int fftSize, Real samplerate, int minimumPeakWidth, int peakCandidateCount)
{
  DEBUG("PITCHACF: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate );

  setFftSize( fftSize, false );
  setSamplerate( samplerate, false );
  setMinimumPeakWidth( minimumPeakWidth, false );
  setPeakCandidateCount( peakCandidateCount, false );
  setup();
}

PitchACF::~PitchACF(){}

void PitchACF::setup(){
  DEBUG("PITCHACF: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;

  _peak.setPeakCount( 1, false );
  _peak.setSortMethod( PeakDetection::BYMAGNITUDE, false );
  _peak.setMinimumPeakWidth( _minimumPeakWidth, false );
  _peak.setCandidateCount( _peakCandidateCount, false );
  _peak.setup();

  _acorr.setInputSize( _halfSize, false );
  _acorr.setMaxLag( _halfSize, false );
  _acorr.setup();

  reset();

  DEBUG("PITCHACF: Finished setup.");
}

void PitchACF::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies){
  _acorr.process(spectrum, &_acorred);
  
  _peak.process(_acorred,
                pitches, saliencies);
  
  _peakInterp.process(_acorred, (*pitches), (*saliencies),
                      pitches, saliencies);

  (*pitches) *= 2.0 * _samplerate / _fftSize;
  
  (*saliencies).cwise() /= _acorred.col(0);
}

void PitchACF::reset(){
  // Initial values

}

int PitchACF::fftSize() const{
  return _fftSize;
}

void PitchACF::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

int PitchACF::minimumPeakWidth() const{
  return _minimumPeakWidth;
}

void PitchACF::setMinimumPeakWidth( int width, bool callSetup ) {
  _minimumPeakWidth = width;
  if ( callSetup ) setup();
}

int PitchACF::peakCandidateCount() const{
  return _peakCandidateCount;
}

void PitchACF::setPeakCandidateCount( int count, bool callSetup ) {
  _peakCandidateCount = count;
  if ( callSetup ) setup();
}

Real PitchACF::samplerate() const{
  return _samplerate;
}
  
void PitchACF::setSamplerate( Real frequency, bool callSetup ){
  _samplerate = frequency;
  if ( callSetup ) setup();
}
