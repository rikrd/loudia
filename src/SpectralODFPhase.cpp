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

#include "SpectralODFPhase.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralODFPhase::SpectralODFPhase(int fftSize, bool weighted, bool normalize) :
  SpectralODFBase(),
  _fftSize( fftSize ),
  _halfSize( _fftSize / 2 + 1 ),
  _weighted( weighted ),
  _normalize( normalize )
{
  
  DEBUG("SPECTRALODFPHASE: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFPhase::~SpectralODFPhase() {}


void SpectralODFPhase::setup() {
  // Prepare the buffers
  DEBUG("SPECTRALODFPHASE: Setting up...");

  _unwrap.setup();

  reset();

  DEBUG("SPECTRALODFPHASE: Finished set up...");
}


void SpectralODFPhase::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("SPECTRALODFPHASE: Processing windowed");
  const int rows = fft.rows();
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);
  
  _unwrap.process(fft.cwise().angle(), &_unwrappedAngle);

  DEBUG("SPECTRALODFPHASE: Processing unwrapped");
  
  phaseDeviation(fft, _unwrappedAngle, odfValue);
  
  DEBUG("SPECTRALODFPHASE: Finished Processing");
}

void SpectralODFPhase::phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  _phaseDiff = spectrumArg.block(1, 0, rows - 1, cols) - spectrumArg.block(0, 0, rows - 1, cols);
  _instFreq = _phaseDiff.block(1, 0, rows - 2, cols) - _phaseDiff.block(0, 0, rows - 2, cols);

  if (_weighted)
    _instFreq.cwise() *= spectrum.block(2, 0, rows - 2, cols).cwise().abs();

  if (_normalize) {
    (*odfValue) = _instFreq.rowwise().sum().cwise() / (cols * spectrum.block(2, 0, rows - 2, cols).cwise().abs().rowwise().sum());
    return;
  }
  
  (*odfValue) = _instFreq.rowwise().sum() / cols;
  return;
}


void SpectralODFPhase::reset() {
  // Initial values
  _unwrap.reset();
}

