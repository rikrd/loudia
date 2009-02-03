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

#include "nmf.h"

#include "utils.h"
#include <Eigen/LU>

using namespace std;

// import most common Eigen types 
using namespace Eigen;

NMF::NMF(int fftSize, int numComponents, int maxIterations, Real maxError, Real eps) :
  _fftSize(fftSize),
  _numComponents(numComponents),
  _maxIterations(maxIterations),
  _maxError(maxError),
  _eps(eps)
{
  
  DEBUG("NMF: Constructor fftSize: " << _fftSize
        << " numComponents: " << _numComponents
        << " maxIterations: " << _maxIterations
        << " maxError: " << _maxError );
  
  setup();
}

NMF::~NMF() {}


void NMF::setup() {
  // Prepare the buffers
  DEBUG("NMF: Setting up...");
  
  reset();

  DEBUG("NMF: Finished set up...");
}


void NMF::process(const MatrixXR& spectrumAbs, MatrixXR* components, MatrixXR* gains) {
  DEBUG("NMF: Processing windowed");
  const int rows = spectrumAbs.rows();
  const int cols = spectrumAbs.cols();
  
  // The X matrix is spectrumAbs.transpose()
  // Some beleive it can be useful to normalize
  
  // The W matrix is (*components).transpose()
  (*components).resize(_numComponents, cols);
  
  // The H matrix is (*gains).transpose()
  (*gains).resize(rows, _numComponents);
  
  DEBUG("NMF: Components resized rows: " << rows << " halfCols: " << halfCols);

  // Initializing W and H
  // TODO: initialize with a Normal distribution
  (*components).setRandom();
  (*components) = (*components).cwise().abs();

  (*gains).setRandom();
  (*gains) = (*gains).cwise().abs();
  
  for (int iter = 0; iter < _maxIterations; iter ++) {
    _xOverWH = spectrumAbs.transpose().cwise() / (((*components).transpose() * (*gains).transpose()).cwise() + _eps);
    
    // Multiplicative update rules of W and H by (Lee and Seung 2001)
    (*components).transpose().cwise() *= (_xOverWH * (*gains)).cwise() / (MatrixXR::Ones(cols, 1) * (*gains).colwise().sum());
    
    (*gains).transpose().cwise() *= ((*components) * _xOverWH).cwise() / ((*components).transpose().colwise().sum().transpose() * MatrixXR::Ones(1, rows));

    // Renormalize so rows of H have constant energy
    _norms = (*gains).colwise().norm();
    
    (*components).transpose().cwise() *= MatrixXR::Ones(cols, 1) * _norms;
    (*gains).transpose().cwise() /= _norms.transpose() * MatrixXR::Ones(1, rows);
  }
  
  DEBUG("NMF: Finished Processing");
}

void NMF::reset() {
  // Initial values
}
