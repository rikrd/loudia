/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

ODFComplex::ODFComplex(int fftLength) : _unwrap((int)(fftLength / 2.0)) {
  
  DEBUG("ODFComplex: Constructor fftLength: " << fftLength);
  
  _fftLength = fftLength;
  
  setup();
}

ODFComplex::~ODFComplex() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void ODFComplex::setup() {
  // Prepare the buffers
  DEBUG("ODFComplex: Setting up...");

  _unwrap.setup();

  reset();

  DEBUG("ODFComplex: Finished set up...");
}


void ODFComplex::process(MatrixXC fft, MatrixXR* odfValue) {
  DEBUG("ODFComplex: Processing windowed");

  (*odfValue).resize(1, 1);
  _spectrum.resize(fft.rows(), (int)ceil(_fftLength / 2.0));
  
  _spectrum = fft.block(0, 0, fft.rows(), (int)ceil(_fftLength / 2.0));

  _unwrap.process(_spectrum.angle().real().cast<Real>(), &_unwrappedAngle);
  
  (*odfValue)(0, 0) = spectralDistanceEuclidean(_spectrum, _spectrum.cwise().abs(), _unwrappedAngle);
  
  DEBUG("ODFComplex: Finished Processing");
}

Real ODFComplex::spectralDistanceEuclidean(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  if (rows < 3) {
    // Throw not enough rows
  }
  
  _spectrumPredict.resize(1, cols);
  _predictionError.resize(1, cols);

  polar(spectrumAbs.row(rows - 2), 2.0*spectrumArg.row(rows - 2) - spectrumArg.row(rows - 3), &_spectrumPredict);
  
  _predictionError = (_spectrumPredict.row(0) - spectrum.row(rows - 1)).cwise().abs();
  
  _predictionError(0,0) = 0.0;
  
  return _predictionError.sum() / (cols-1) * sqrt(2.0);
}

Real ODFComplex::spectralDistanceEuclideanWeighted(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg) {
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

Real ODFComplex::spectralDistanceHypot(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg) {
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
