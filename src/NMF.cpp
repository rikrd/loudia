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

#include "NMF.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

NMF::NMF(int vectorSize, int factor, int maxIterations, Real eps) :
  _vectorSize(vectorSize),
  _factor(factor),
  _maxIterations(maxIterations),
  _eps(eps)
{
  
  DEBUG("NMF: Constructor vectorSize: " << _vectorSize
        << " factor: " << _factor
        << " maxIterations: " << _maxIterations );
  
  setup();
}

NMF::~NMF() {}


void NMF::setup() {
  // Prepare the buffers
  DEBUG("NMF: Setting up...");
  
  reset();

  DEBUG("NMF: Finished set up...");
}


void NMF::process(const MatrixXR& v, MatrixXR* w, MatrixXR* h) {
  DEBUG("NMF: Processing ...");
  const int rows = v.rows();
  const int cols = v.cols();
  
  // The X matrix is v.transpose()
  // Some beleive it can be useful to normalize
  
  // The W matrix is (*w).transpose()
  (*w).resize(_factor, cols);
  
  // The H matrix is (*h).transpose()
  (*h).resize(rows, _factor);
  
  // Initializing W and H
  // TODO: initialize with a Normal distribution
  (*w).setRandom();
  (*w) = (*w).cwise().abs();

  (*h).setRandom();
  (*h) = (*h).cwise().abs();
  
  for (int iter = 0; iter < _maxIterations; iter ++) {
    _xOverWH = v.transpose().cwise() / (((*w).transpose() * (*h).transpose()).cwise() + _eps);
    
    // Multiplicative update rules of W and H by (Lee and Seung 2001)
    (*w).transpose().cwise() *= (_xOverWH * (*h)).cwise() / (MatrixXR::Ones(cols, 1) * (*h).colwise().sum());
    
    (*h).transpose().cwise() *= ((*w) * _xOverWH).cwise() / ((*w).transpose().colwise().sum().transpose() * MatrixXR::Ones(1, rows));

    // Renormalize so rows of H have constant energy
    _norms = (*h).colwise().norm();
    
    (*w).transpose().cwise() *= MatrixXR::Ones(cols, 1) * _norms;
    (*h).transpose().cwise() /= _norms.transpose() * MatrixXR::Ones(1, rows);
  }
  
  DEBUG("NMF: Finished Processing");
}

void NMF::reset() {
  // Initial values
}
