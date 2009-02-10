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
#include <vector>

#include "chebyshev.h"
#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

Chebyshev::Chebyshev( int order, Real freq, Real rippleDB, int channels, ChebyshevType chebyshevType) : 
  _order(order),
  _freq(freq),
  _rippleDB(rippleDB),
  _channels(channels),
  _filter(channels),
  _chebyshevType(chebyshevType)
                                                                                           
{
  DEBUG("CHEBYSHEV: Constructor order: " << order << 
        ", freq: " << freq << 
        ", rippleDB: " << rippleDB );

  if ( order < 1 ) {
    // Throw an exception
  }
  
  setup();
  
  DEBUG("CHEBYSHEV: Constructed");
}

void Chebyshev::chebyshev1(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  (*zeros) = MatrixXC::Zero(channels, 1);
 
  Real eps = sqrt(pow(10, (0.1 * rippleDB)) - 1.0);

  MatrixXC n;
  range(1, order + 1, order, channels, &n);

  Real mu = 1.0 / order * log((1.0 + sqrt( 1 + eps * eps)) / eps);

  MatrixXC theta = ((n * 2).cwise() - 1.0) / order * M_PI / 2.0;

  (*poles) = -sinh(mu) * theta.cwise().sin() + Complex(0, 1) * cosh(mu) * theta.cwise().cos();

  Complex gainComplex = 1.0;

  for ( int i = 0; i < (*poles).cols(); i++ ) {
    gainComplex *= -(*poles)(0, i);
  }

  (*gain) = gainComplex.real();
  
  if ( order % 2 == 0 ) {
    (*gain) /= sqrt((1 + eps * eps));
  }
}

void Chebyshev::chebyshev2(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  Real de = 1.0 / sqrt(pow(10, (0.1 * rippleDB)) - 1.0);
  Real mu = arcsinh(1.0 / de) / order;
  
  int m;
  MatrixXC n(channels, order);
  if(order % 2) {
    m = order - 1;
    MatrixXC nFirst;
    range(1, order - 1, order/2 - 1, channels, &nFirst);

    MatrixXC nSecond;
    range(order + 2, 2 * order, order/2 - 1, channels, &nSecond);
    
    n << nFirst, nSecond;
  } else{
    m = order;
    range(1, 2*order, order, channels, &n);
  }
  
  (*zeros) = (Complex(0,1) / ((n + M_PI) / (2.0*order)).cwise().cos()).conj();

  MatrixXC rng;
  range(1, 2*order, order, channels, &rng);
  
  (*poles) = (Complex(0,1) * (M_PI * rng / (2.0*order) + M_PI/2.0)).cwise().exp();
  (*poles) = (sinh( mu ) * (*poles).real() + Complex(0, 1) * cosh( mu ) * (*poles).imag()).cwise().inverse();
  
  // TODO: use the Eigen prod() when the patch is accepted
  Real num = 1.0;
  //k = (numpy.prod(-p,axis=0)/numpy.prod(-z,axis=0)).real
}


void Chebyshev::setup(){
  DEBUG("CHEBYSHEV: Setting up...");

  DEBUG("CHEBYSHEV: Getting zpk");  
  // Get the chebyshev z, p, k
  MatrixXC zeros, poles;
  Real gain;

  switch( _chebyshevType ){
  case I:
    chebyshev1(_order, _rippleDB, _channels, &zeros, &poles, &gain);
    break;

  case II:
    chebyshev2(_order, _rippleDB, _channels, &zeros, &poles, &gain);
    break;
  }
  
  DEBUG("CHEBYSHEV: zeros:" << zeros );
  DEBUG("CHEBYSHEV: poles:" << poles );
  DEBUG("CHEBYSHEV: gain:" << gain );

  // Convert zpk to ab coeffs
  MatrixXC a;
  MatrixXC b;
  zpkToCoeffs(zeros, poles, gain, &b, &a);
  
  // Since we cannot create matrices of Nx0
  // we have created at least one Zero in 0
  // Now we must remove the last coefficient from b

  MatrixXC temp = b.block(0, 0, b.rows(), b.cols()-1);
  b = temp;

  // Get the warped critical frequency
  Real fs = 2.0;
  Real warped = 2.0 * fs * tan( M_PI * _freq / fs );

  // Warpped coeffs
  MatrixXC wa;
  MatrixXC wb;
  lowPassToLowPass(b, a, warped, &wb, &wa);
 
  // Digital coeffs
  MatrixXR da;
  MatrixXR db;
  bilinear(wb, wa, fs, &db, &da);
  
  // Set the coefficients to the filter
  _filter.setB( db.transpose() );
  _filter.setA( da.transpose() );
  _filter.setup();
  
  DEBUG("CHEBYSHEV: Finished set up...");
}

void Chebyshev::a(MatrixXR* a) {
  _filter.a(a);
}

void Chebyshev::b(MatrixXR* b) {
  _filter.b(b);
}

void Chebyshev::process(MatrixXR samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void Chebyshev::reset(){
  // Initial values
  _filter.reset();
}
