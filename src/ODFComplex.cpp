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

#include "ODFComplex.h"
#include "Unwrap.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

ODFComplex::ODFComplex(int fftSize, bool rectified) :
  ODFBase(),
  _fftSize( fftSize ),
  _halfSize( _fftSize / 2 + 1 ),
  _rectified( rectified ),
  _unwrap( _halfSize )
{
  
  DEBUG("ODFCOMPLEX: Constructor fftSize: " << _fftSize);
  
  setup();
}

ODFComplex::~ODFComplex() {}


void ODFComplex::setup() {
  // Prepare the buffers
  DEBUG("ODFCOMPLEX: Setting up...");

  _unwrap.setup();

  reset();

  DEBUG("ODFCOMPLEX: Finished set up...");
}


void ODFComplex::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFCOMPLEX: Processing windowed");
  const int rows = fft.rows();
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);

  _unwrap.process(fft.cwise().angle(), &_unwrappedAngle);

  DEBUG("ODFCOMPLEX: Processing unwrapped");
  
  spectralDistanceEuclidean(fft, fft.cwise().abs(), _unwrappedAngle, odfValue);
  
  DEBUG("ODFCOMPLEX: Finished Processing");
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
