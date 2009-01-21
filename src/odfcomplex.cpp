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

#include "odfcomplex.h"
#include "unwrap.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFComplex::ODFComplex(int fftLength, bool rectified) :
  ODFBase(),
  _fftLength(fftLength),
  _rectified(rectified),
  _unwrap((int)(fftLength / 2.0))
{
  
  DEBUG("ODFComplex: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFComplex::~ODFComplex() {}


void ODFComplex::setup() {
  // Prepare the buffers
  DEBUG("ODFComplex: Setting up...");

  _unwrap.setup();

  reset();

  DEBUG("ODFComplex: Finished set up...");
}


void ODFComplex::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFComplex: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);
  _spectrum.resize(rows, halfCols);

  DEBUG("ODFComplex: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrum = fft.block(0, 0, rows, halfCols);

  DEBUG("ODFComplex: Specturum halved");

  _unwrap.process(_spectrum.cwise().angle(), &_unwrappedAngle);

  DEBUG("ODFComplex: Processing unwrapped");
  
  spectralDistanceEuclidean(_spectrum, _spectrum.cwise().abs(), _unwrappedAngle, odfValue);
  
  DEBUG("ODFComplex: Finished Processing");
}

void ODFComplex::spectralDistanceEuclidean(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg, MatrixXR* odfValue) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  _spectrumPredict.resize(rows - 2, cols);
  _predictionError.resize(rows - 2, cols);

  polar(spectrumAbs.block(1, 0, rows - 2, cols), 2.0 * spectrumArg.block(1, 0, rows - 2, cols) - spectrumArg.block(0, 0, rows - 2, cols), &_spectrumPredict);
  
  _predictionError = (_spectrumPredict - spectrum.block(0, 0, rows - 2, cols)).cwise().abs();

  if (_rectified) {
    _predictionError = (_spectrumPredict.cwise().abs().cwise() <= spectrum.block(0, 0, rows - 2, cols).cwise().abs()).select(_predictionError, 0.0);
  }
  
  //_predictionError.col(0) = 0.0;
  
  (*odfValue) = _predictionError.rowwise().sum() / cols;
  return;
}

void ODFComplex::reset() {
  // Initial values
  _unwrap.reset();
}
