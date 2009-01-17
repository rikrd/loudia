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

#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/QR> 

#include <cmath>

#include "utils.h"

using namespace std;

void roots(const MatrixXR& poly, MatrixXC* result) {
  const int coeffs = poly.cols();
  
  if ( coeffs <= 1 ) {
    // Throw error about wrong input length
  }
  
  // Prepare the output
  (*result).resize(coeffs - 1, 1);

  // Build companion matrix and find its eigenvalues (the roots)
  MatrixXR companion(coeffs - 1, coeffs - 1);
  companion.corner( Eigen::TopRight, coeffs - 2, 1 ).setZero();
  companion.corner( Eigen::BottomLeft, coeffs - 2, coeffs - 2).setIdentity();
  companion.row(0) = -poly.corner( Eigen::TopRight, 1, coeffs - 1 ) / poly(0, 0);
  
  // Get the eigen values
  (*result) = Eigen::EigenSolver<MatrixXR>(companion).eigenvalues();
}

void poly(const MatrixXR&  roots, MatrixXC* result) {
  const int nroots = roots.cols();
  
  // Prepare the output
  (*result).resize(1, nroots + 1);
}

void reverseCols(MatrixXC* in) { 
  const int cols = (*in).cols();
  
  for(int i = 0; i < cols / 2; i++ ){
    (*in).col(i).swap((*in).col(cols - i - 1));
  }
}

void reverseCols(MatrixXR* in) {
  const int cols = (*in).cols();
  
  for(int i = 0; i < cols / 2; i++ ){
    (*in).col(i).swap((*in).col(cols - i - 1));
  }
}

void rowCumsum(MatrixXR* in) { 
  const int rows = (*in).rows();
  
  for(int i = 1; i < rows; i++ ){
    (*in).row(i) += (*in).row(i-1);
  }
}


void polar(const MatrixXR&  mag, const MatrixXR&  phase, MatrixXC* complex) {
  if ((mag.rows() != phase.rows()) || (mag.cols() != phase.cols())) {
    // Throw an error
  }

  (*complex).resize(mag.rows(), mag.cols());

  for(int i = 0; i < mag.rows(); i++){
    for(int j = 0; j < mag.cols(); j++){
      (*complex)(i, j) = polar(mag(i, j), phase(i, j));
    }
  }
}

void coeffsToZpk(const MatrixXR&  b, const MatrixXR&  a, MatrixXC* zeros, MatrixXC* poles, Real* gain){
  (*gain) = b(0, 0);
  MatrixXR bTemp = b;
  bTemp /= b(0, 0);
  roots(bTemp, zeros);
  roots(a, poles);
}

Real asinc(int M, Real omega) {
  return omega == 0 ? M : sin(M * omega / 2.0) / sin(omega / 2.0) ;
}


void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        int mainLobeBandwith, Real alpha, Real beta, MatrixXR* spectrum) {
  
  const int begin = max((position - mainLobeBandwith / 2.0), 0.0);
  const int end = min(ceil(position + mainLobeBandwith / 2.0 + 1), fftSize/2.0);

  spectrum->resize(1, end - begin);
  
  const Real omegaM = 2.0 * M_PI / fftSize;
  
  for ( int row = 0; row < spectrum->rows(); row++ ) { 
    for ( int i = begin; i < end; i++ ) {
      
      Real omega = 2.0 * M_PI * (i - position) / fftSize;

      (*spectrum)(row, i-begin) = magnitude * abs(alpha * asinc(windowSize, omega) + beta * (asinc(windowSize, omega - omegaM) + asinc(windowSize, omega + omegaM)));

    }
  }
}

void hannTransform(Real position, Real magnitude, 
                   int windowSize, int fftSize,
                   int mainLobeBandwith, MatrixXR* spectrum) {

  raisedCosTransform(position, magnitude, windowSize, fftSize, mainLobeBandwith, 0.5, 0.5, spectrum);

}


void hammingTransform(Real position, Real magnitude, 
                      int windowSize, int fftSize,
                      int mainLobeBandwith, MatrixXR* spectrum) {

  raisedCosTransform(position, magnitude, windowSize, fftSize, mainLobeBandwith, 0.53836, 0.46164, spectrum);

}



/*
void zpkToCoeffs(MatrixXC zeros, MatrixXC poles, Real gain, MatrixXC* b, MatrixXC* a):
    """Return polynomial transfer function representation from zeros
    and poles

    Inputs:

      z, p --- sequences representing the zeros and poles.
      k --- system gain.

    Outputs: (b,a)

      b, a --- numerator and denominator polynomials.
    """
    z = atleast_1d(z)
    k = atleast_1d(k)
    if len(z.shape) > 1:
        temp = poly(z[0])
        b = zeros((z.shape[0], z.shape[1]+1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]]*z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * poly(z[i])
    else:
        b = k * poly(z)
    a = poly(p)
    return b, a
*/
