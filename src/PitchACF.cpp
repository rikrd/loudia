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

#include "PitchACF.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

PitchACF::PitchACF(int fftSize, Real samplerate, int minPeakWidth, int peakCandidateCount) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _samplerate( samplerate ),
  _peak(1, PeakDetect::BYMAGNITUDE, minPeakWidth, peakCandidateCount),
  _peakInterp(),
  _acorr(_halfSize, _halfSize)
{
  DEBUG("PITCHACF: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate );

  setup();
}

PitchACF::~PitchACF(){}

void PitchACF::setup(){
  DEBUG("PITCHACF: Setting up...");

  reset();

  DEBUG("PITCHACF: Finished setup.");
}

void PitchACF::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies){
  _acorr.process(spectrum, &_acorred); 
 
  _acorredC = _acorred.cast<Complex>();
  
  _peak.process(_acorredC,
                pitches, saliencies, &_phases);
  
  _peakInterp.process(_acorredC, (*pitches), (*saliencies), _phases,
                      pitches, saliencies, &_phases);

  (*pitches) *= _samplerate / _fftSize;

  (*saliencies).cwise() /= _acorred.col(0);
}

void PitchACF::reset(){
  // Initial values

}
