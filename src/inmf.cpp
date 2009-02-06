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

#include "inmf.h"

#include "utils.h"
#include <Eigen/LU>

using namespace std;

// import most common Eigen types 
using namespace Eigen;

INMF::INMF(int fftSize, int numComponents, int numPast, Real pastCoeff, Real newCoeff, int maxIterations, Real maxError, Real eps) :
  _fftSize(fftSize),
  _numComponents(numComponents),
  _maxIterations(maxIterations),
  _maxError(maxError),
  _eps(eps),
  _numPast(numPast),
  _pastCoeff(pastCoeff),
  _newCoeff(newCoeff)
{
  
  DEBUG("INMF: Constructor fftSize: " << _fftSize
        << " numComponents: " << _numComponents
        << " numPast: " << _numPast
        << " pastCoeff: " << _pastCoeff
        << " newCoeff: " << _newCoeff
        << " maxIterations: " << _maxIterations
        << " maxError: " << _maxError );
  
  setup();
}

INMF::~INMF() {}


void INMF::setup() {
  // Prepare the buffers
  DEBUG("INMF: Setting up...");
  
  _V.resize(_numPast, _fftSize);
  _H.resize(_numPast, _numComponents);
  _W.resize(_numComponents, _fftSize);
  
  reset();

  DEBUG("INMF: Finished set up...");
}


void INMF::process(const MatrixXR& v, MatrixXR* w, MatrixXR* h) {
  DEBUG("INMF: Processing windowed");
  const int rows = v.rows();
  const int cols = v.cols();
  
  if ( cols != _numComponents ) {
    // Throw error wrong number of columns
  }

  // The X matrix is spectrumAbs.transpose()
  // Some beleive it can be useful to normalize
  
  // The W matrix is (*w).transpose()
  (*w) = _W;
  
  // The H matrix is (*h).transpose()
  (*h).resize(rows, _numComponents);
  
  // Initializing h
  // TODO: initialize with a Normal distribution
  (*h).setRandom();
  (*h) = (*h).cwise().abs();

  DEBUG("INMF: Begin rows");
  for (int row = 0; row < rows; row++ ) {

    // Calculate beta * VHt
    _VH = _pastCoeff * (_V.transpose() * _H);

    // Calculate beta * HHt
    _HH = _pastCoeff * (_H.transpose() * _H);

    DEBUG("INMF: Begin iterations");
    for ( int iter = 0; iter < _maxIterations; iter++ ) {
      /*
      MatrixXR Wv = (*w) * v.row(row).transpose();
      MatrixXR WWh = ((*w) * (*w).transpose()) * (*h).row(row).transpose();
      for ( int a = 0; a < _numComponents; a++ ) {
        (*h).row( row )(a) *= Wv(a) / WWh(a);
      }
      */

      //DEBUG("INMF: Update h");
      // Eq. 9 in Bucak 2008
      /*
      cout << "---------------" << endl;
      cout << "(*w): " << (*w).rows() << " , " << (*w).cols() << endl;
      cout << "---------------" << endl;
      cout << "(*h): " << (*h).rows() << " , " << (*h).cols() << endl;
      cout << "---------------" << endl;
      cout << "v: " << v.rows() << " , " << v.cols() << endl;
      cout << "---------------" << endl;
      */
      (*h).row( row ).cwise() *= (((*w) * v.row(row).transpose()).cwise() / ((((*w) * (*w).transpose()) * (*h).row( row ).transpose()).cwise() + _eps)).transpose();

      //DEBUG("INMF: Update W");
      // Eq. 13 in Bucak 2008
      /*
      cout << "---------------" << endl;
      cout << "(*w): " << (*w).rows() << " , " << (*w).cols() << endl;
      cout << "_HH: " << _HH.rows() << " , " << _HH.cols() << endl;
      cout << "(*h): " << (*h).rows() << " , " << (*h).cols() << endl;
      cout << "_VH: " << _VH.rows() << " , " << _VH.cols() << endl;
      cout << "---------------" << endl;

      cout << ((*w).transpose() * _HH) << endl;
      cout << "============" << endl;
      cout << (_newCoeff * ((*h).row( row ).transpose() * (*h).row( row ))) << endl;
      cout << "**********" << endl;
      */
      (*w).cwise() *= ((_VH + (_newCoeff * v.row( row ).transpose() * (*h).row( row ))).cwise() / (((*w).transpose() * (_HH + _newCoeff * (*h).row( row ).transpose() * (*h).row( row ))).cwise() + _eps)).transpose();
    }

    DEBUG("INMF: Shift and update H");
    // Update the past H
    // TODO: when Eigen will have rotate use this
    //_H.rowwise().rotate(-1);
    rowShift(&_H, -1);
    _H.row(_numPast - 1) = (*h).row( row );

    DEBUG("INMF: Shift and update V");
    // Update the past V
    // TODO: when Eigen will have rotate use this
    //_V.rowwise().rotate(-1);
    rowShift(&_V, -1);
    _V.row(_numPast - 1) = v.row( row );
    
    DEBUG("INMF: Keep W");
    // Keep the past W
    _W = (*w);
  }

  DEBUG("INMF: Finished Processing");
}

void INMF::reset() {
  // Initial W, H and V
  // TODO: initialize with a Normal distribution
  _W.setRandom();
  _W = _W.cwise().abs();

  _H.setRandom();
  _H = _H.cwise().abs();

  _V.setZero();
  
}
