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

#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/QR> 

#include <iostream>
#include <cmath>
#include <vector>

#include "chebyshev.h"

#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

void Chebyshev::reverse(MatrixXC* in){
  const int N = (*in).cols();
  
  for(int i = 0; i < N / 2; i++ ){
    (*in).col(i).swap((*in).col(N - i - 1));

  }
}

void Chebyshev::reverse(MatrixXR* in){
  const int N = (*in).cols();
  
  for(int i = 0; i < N / 2; i++ ){
    (*in).col(i).swap((*in).col(N - i - 1));

  }
}


void Chebyshev::roots(MatrixXR poly, MatrixXC* roots){
  const int N = poly.cols();
  (*roots).resize(1, N-1);
  
  if ( N > 1 ) {
    // Build companion matrix and find its eigenvalues (the roots)
    MatrixXR A = MatrixXR::Zero(N - 1, N - 1);
    A.corner( Eigen::BottomLeft, N - 2, N - 2).diagonal().setOnes();
    A.row(0) = -poly.corner( Eigen::TopRight, 1, N - 1 ) / poly(0, 0);
    
    // Get the eigen values
    (*roots).set(Eigen::EigenSolver<MatrixXR>(A).eigenvalues().transpose());
  }
  
  reverse(roots);
  
}

/*
void Chebyshev::tf2zpk(MatrixXC b, MatrixXC a, MatrixXC* zeros, MatrixXC* poles, Real* gain){
  // Return zero, pole, gain (z,p,k) representation from a numerator,
  // denominator representation of a linear filter.
  (*gain) = b(0);
  b /= b[0];
  z = roots(b);
  p = roots(a);
}

def zpk2tf(MatrixXC zeros, MatrixXC poles, Real gain, MatrixXC* b, MatrixXC* a):
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

Chebyshev::Chebyshev(int channels, int order, Real rippleDB, Real samplerate) : _filter(channels)
{
  DEBUG("CHEBYSHEV: Constructor order: " << order << ", rippleDB: " << rippleDB << ", samplerate: " << samplerate);

  if ( order < 1 ) {
    // Throw an exception
  }

  if ( samplerate <= 0 ) {
    // Throw an exception
  }
  
  _channels = channels;
  _order = order;
  _rippleDB = rippleDB;
  _samplerate = samplerate;

  setup();
  DEBUG("CHEBYSHEV: Constructed");
}

void Chebyshev::setup(){
  DEBUG("CHEBYSHEV: Setting up...");
  
  MatrixXC zeros(1,1);
  zeros.setZero();
  
  Real eps = sqrt(pow(10, (0.1 * _rippleDB) - 1.0));

  MatrixXR n(1, _order + 1);
  for ( int i = 0; i < n.cols(); i++ ) {
    n(0, i) = i;
  }
  
  Real mu = 1.0 / _order * log((1.0 + sqrt( 1 + eps * eps)) / eps);

  MatrixXC theta = ((n * 2).cwise() - 1.0) / _order * M_PI / 2.0;

  MatrixXC poles = -sinh(mu) * theta.cwise().sin() + Complex(0, 1) * cosh(mu) * theta.cwise().cos();

  Complex gainComplex = 1.0;
  for ( int i = 0; i < poles.cols(); i++ ) {
    gainComplex *= -poles(0, i);
  }

  Real gain = gainComplex.real();
  
  if ( _order % 2 == 0 ) {
    gain = gain / sqrt((1 + eps * eps));
  }
  
  DEBUG("CHEBYSHEV: zeros:" << zeros );
  DEBUG("CHEBYSHEV: poles:" << poles );
  DEBUG("CHEBYSHEV: gain:" << gain );

  DEBUG("CHEBYSHEV: Finished set up...");
}

void Chebyshev::process(MatrixXR samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void Chebyshev::reset(){
  // Initial values
  _filter.reset();
}
