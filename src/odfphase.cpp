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

#include "typedefs.h"
#include "debug.h"

#include <cmath>

#include "odfphase.h"
#include "unwrap.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFPhase::ODFPhase(int fftLength, bool weighted, bool normalize) :
  ODFBase(),
  _fftLength(fftLength),
  _weighted(weighted),
  _normalize(normalize),
  _unwrap((int)(fftLength / 2.0))
{
  
  DEBUG("ODFPhase: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFPhase::~ODFPhase() {}


void ODFPhase::setup() {
  // Prepare the buffers
  DEBUG("ODFPhase: Setting up...");

  _unwrap.setup();

  reset();

  DEBUG("ODFPhase: Finished set up...");
}


void ODFPhase::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFPhase: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);

  DEBUG("ODFPhase: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrum = fft.block(0, 0, rows, halfCols);

  DEBUG("ODFPhase: Specturum halved");

  _unwrap.process(_spectrum.cwise().angle(), &_unwrappedAngle);

  DEBUG("ODFPhase: Processing unwrapped");
  
  phaseDeviation(_spectrum, _unwrappedAngle, odfValue);
  
  DEBUG("ODFPhase: Finished Processing");
}

void ODFPhase::phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue) {
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


void ODFPhase::reset() {
  // Initial values
  _unwrap.reset();
}
