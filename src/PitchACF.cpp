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

PitchACF::PitchACF(int fftSize, Real samplerate) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _samplerate( samplerate ),
  _peak(1, false, 40),
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
  const int rows = spectrum.rows();

  (*pitches).resize( rows, 1 );
  (*saliencies).resize( rows, 1 );

  MatrixXR phases;

  _acorr.process(spectrum, &_acorred); 
 
  _acorredC = _acorred.cast<Complex>();

  _peak.process(_acorredC,
                pitches, saliencies, &phases);
  
  _peakInterp.process(_acorredC, (*pitches), (*saliencies), phases,
                      pitches, saliencies, &phases);

  (*pitches) *= _samplerate / _fftSize;
}

void PitchACF::reset(){
  // Initial values

}
