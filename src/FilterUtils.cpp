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

void chebyshev1(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
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

void chebyshev2(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain) {
  Real de = 1.0 / sqrt(pow(10, (0.1 * rippleDB)) - 1.0);
  Real mu = asinh(1.0 / de) / order;

  MatrixXC n;
  if(order % 2) {
    n.resize(channels, order - 1);
    MatrixXC nFirst;
    range(1, order , order/2, channels, &nFirst);

    MatrixXC nSecond;
    range(order + 2, (2 * order) + 1, order/2, channels, &nSecond);
    
    n << nFirst, nSecond;
  } else{
    n.resize(channels, order);
    range(1, (2 * order) + 1, order, channels, &n);
  }
  
  (*zeros) = (Complex(0,1) * ((n * M_PI) / (2.0 * order)).cwise().cos().cwise().inverse()).conjugate();

  MatrixXC rng;
  range(1, (2 * order) + 1, order, channels, &rng);
  
  (*poles) = (Complex(0,1) * (((M_PI * rng) / (2.0*order)).cwise() + M_PI / 2.0)).cwise().exp();

  (*poles) = (((*poles).real().cast<Complex>() * sinh( mu )) + (Complex(0, 1) * cosh( mu ) * (*poles).imag().cast<Complex>())).cwise().inverse();

  // TODO: gain should be a vector (one gain per channel)
  (*gain) = ((-(*poles)).rowwise().prod().cwise() / (-(*zeros)).rowwise().prod()).real().sum();
}
