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

#include "NMF.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

NMF::NMF(int inputSize, int componentCount, int iterationCount, Real epsilon)
{
  
  LOUDIA_DEBUG("NMF: Constructor inputSize: " << inputSize
        << " componentCount: " << componentCount
        << " iterationCount: " << iterationCount );
  
  setInputSize( inputSize, false );
  setComponentCount( componentCount, false );
  setIterationCount( iterationCount, false );
  setEpsilon( epsilon, false );

  setup();
}

NMF::~NMF() {}


void NMF::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("NMF: Setting up...");
  
  reset();

  LOUDIA_DEBUG("NMF: Finished set up...");
}


void NMF::process(const MatrixXR& v, MatrixXR* w, MatrixXR* h) {
  LOUDIA_DEBUG("NMF: Processing ...");
  const int rows = v.rows();
  const int cols = v.cols();
  
  // The X matrix is v.transpose()
  // Some beleive it can be useful to normalize
  
  // The W matrix is (*w).transpose()
  (*w).resize(_componentCount, cols);
  
  // The H matrix is (*h).transpose()
  (*h).resize(rows, _componentCount);
  
  // Initializing W and H
  // TODO: initialize with a Normal distribution
  (*w).setRandom();
  (*w) = (*w).cwise().abs();

  (*h).setRandom();
  (*h) = (*h).cwise().abs();
  
  for (int iter = 0; iter < _iterationCount; iter ++) {
    _xOverWH = v.transpose().cwise() / (((*w).transpose() * (*h).transpose()).cwise() + _epsilon );
    
    // Multiplicative update rules of W and H by (Lee and Seung 2001)
    (*w).transpose().cwise() *= (_xOverWH * (*h)).cwise() / (MatrixXR::Ones(cols, 1) * (*h).colwise().sum());
    
    (*h).transpose().cwise() *= ((*w) * _xOverWH).cwise() / ((*w).transpose().colwise().sum().transpose() * MatrixXR::Ones(1, rows));

    // Renormalize so rows of H have constant energy
    _norms = (*h).colwise().norm();
    
    (*w).transpose().cwise() *= MatrixXR::Ones(cols, 1) * _norms;
    (*h).transpose().cwise() /= _norms.transpose() * MatrixXR::Ones(1, rows);
  }
  
  LOUDIA_DEBUG("NMF: Finished Processing");
}

void NMF::reset() {
  // Initial values
}

int NMF::inputSize() const {
  return _inputSize;
}
  
void NMF::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}

int NMF::componentCount() const {
  return _componentCount;
}
  
void NMF::setComponentCount( int count, bool callSetup ) {
  _componentCount = count;
  if ( callSetup ) setup();
}

int NMF::iterationCount() const {
  return _iterationCount;
}
  
void NMF::setIterationCount( int count, bool callSetup ) {
  _iterationCount = count;
  if ( callSetup ) setup();
}

Real NMF::epsilon() const {
  return _epsilon;
}
  
void NMF::setEpsilon( Real epsilon, bool callSetup ) {
  _epsilon = epsilon;
  if ( callSetup ) setup();
}
