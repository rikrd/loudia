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

ODFComplex::ODFComplex(int fftLength) :
  _unwrap((int)(fftLength / 2.0)),
  _fftLength(fftLength),
  ODFBase()
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

  (*odfValue).resize(1, 1);
  _spectrum.resize(fft.rows(), min((int)ceil(_fftLength / 2.0), fft.cols()));

  DEBUG("ODFComplex: Spectrum resized fft.rows(): " << fft.rows() << " (int)ceil(_fftLength / 2.0): " << (int)ceil(_fftLength / 2.0));
  
  _spectrum = fft.block(0, 0, fft.rows(), min((int)ceil(_fftLength / 2.0), fft.cols()));

  DEBUG("ODFComplex: Specturum halved");

  _unwrap.process(_spectrum.cwise().angle().real().cast<Real>(), &_unwrappedAngle);

  DEBUG("ODFComplex: Processing unwrapped");
  
  (*odfValue)(0, 0) = spectralDistanceEuclidean(_spectrum, _spectrum.cwise().abs(), _unwrappedAngle);
  
  DEBUG("ODFComplex: Finished Processing");
}

Real ODFComplex::spectralDistanceEuclidean(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  if (rows < 3) {
    // Throw ValueError not enough rows
  }
  
  _spectrumPredict.resize(1, cols);
  _predictionError.resize(1, cols);

  polar(spectrumAbs.row(rows - 2), 2.0*spectrumArg.row(rows - 2) - spectrumArg.row(rows - 3), &_spectrumPredict);
  
  _predictionError = (_spectrumPredict.row(0) - spectrum.row(rows - 1)).cwise().abs();
  
  _predictionError(0,0) = 0.0;
  
  return _predictionError.sum() / (cols-1) * sqrt(2.0);
}

Real ODFComplex::spectralDistanceEuclideanWeighted(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  if (rows < 3) {
    // Throw not enough rows
  }

  _spectrumPredict.resize(1, cols);
  _predictionError.resize(1, cols);

  polar(spectrumAbs.row(rows - 2), 2.0*spectrumArg.row(rows - 2) - spectrumArg.row(rows - 3), &_spectrumPredict);

  _predictionError = ((_spectrumPredict.row(0) - spectrum.row(rows - 1)) * spectrumAbs.row(rows - 1)).cwise().abs() / spectrumAbs.row(rows - 1).sum();

  _predictionError(0,0) = 0.0;
  
  return _predictionError.sum() / (cols-1);
}

Real ODFComplex::spectralDistanceHypot(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  if (rows < 3) {
    // Throw not enough rows
  }
  
  _spectrumPredict.resize(1, cols);
  _predictionError.resize(1, cols);

  polar(spectrumAbs.row(rows - 2), 2.0*spectrumArg.row(rows - 2) - spectrumArg.row(rows - 3), &_spectrumPredict);
  
  MatrixXR _predictionErrorReal = (_spectrumPredict.row(0) - spectrum.row(rows - 1)).real();
  MatrixXR _predictionErrorImag = (_spectrumPredict.row(0) - spectrum.row(rows - 1)).imag();

  
  _predictionError(0,0) = 0.0;
  
  return _predictionError.sum() / (cols-1);
}


void ODFComplex::reset() {
  // Initial values
  _unwrap.reset();
}
