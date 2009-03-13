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

#include "FilterUtils.h"

using namespace std;

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


void coeffsToZpk(const MatrixXR&  b, const MatrixXR&  a, MatrixXC* zeros, MatrixXC* poles, Real* gain){
  (*gain) = b(0, 0);
  MatrixXR bTemp = b;
  bTemp /= b(0, 0);
  roots(bTemp, zeros);
  roots(a, poles);
}

void zpkToCoeffs(const MatrixXC& zeros, const MatrixXC& poles, Real gain, MatrixXC*  b, MatrixXC*  a){
  poly( zeros, b );
  (*b) *= gain;
  
  poly( poles, a );
}

void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();

  const int rows = a.rows();

  const int maxsize = max(asize, bsize);    

  *aout = a;
  *bout = b;

  MatrixXR pwo;
  range(maxsize-1, -1, maxsize, rows, &pwo);

  pwo = pwo.cwise().expN( freq );
 
  int start1 = max(bsize - asize, 0);
  int start2 = max(asize - bsize, 0);
  
  for ( int i = 0; i < bsize; i++ ) {
    (*bout).col(i).cwise() *= pwo.col( start2 + i ).cwise().inverse() * pwo.col( start1 );
  }

  for ( int i = 0; i < asize; i++ ) {
    (*aout).col(i).cwise() *= pwo.col( start1 + i ).cwise().inverse() * pwo.col( start1 );
  }

  normalize((*bout), (*aout));
}


void lowPassToHighPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();

  const int rows = a.rows();

  const int maxsize = max(asize, bsize);    

  (*aout) = MatrixXC::Zero(rows, maxsize);
  (*bout) = MatrixXC::Zero(rows, maxsize);

  (*aout).block(0, 0, rows, asize) = a.rowwise().reverse();
  (*bout).block(0, 0, rows, bsize) = b.rowwise().reverse();

  MatrixXR pwo;
  range(0, maxsize, maxsize, rows, &pwo);

  pwo = pwo.cwise().expN( freq );
  
  (*aout) = (*aout).cwise() * pwo;
  (*bout) = (*bout).cwise() * pwo;

  normalize((*bout), (*aout));
}


void lowPassToBandPass(const MatrixXC& b, const MatrixXC& a, Real freq, Real bandwidth, MatrixXC*  bout, MatrixXC*  aout) { 
  const int asize = a.cols();
  const int bsize = b.cols();

  const int rows = a.rows();

  const int maxsize = max(asize - 1, bsize - 1);
  
  const int maxasize = maxsize + asize - 1;
  const int maxbsize = maxsize + bsize - 1;

  const Real freqSq = freq * freq;

  (*aout) = MatrixXC::Zero(rows, maxasize + 1);
  (*bout) = MatrixXC::Zero(rows, maxbsize + 1);

  for ( int j = 0; j < (maxbsize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < (i + 1); k++ ) {
        if ( maxsize - i + 2 * k == j ) {
          val += (Real)comb(i, k) * b.col(bsize - 1 - i) * pow(freqSq, (Real)(i-k)) / pow(bandwidth, (Real)i);
        }
      }
    }
    
    (*bout).col(maxbsize - j) = val;
  }

  for ( int j = 0; j < (maxasize + 1); j++ ) {
    MatrixXC val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < (i + 1); k++ ) {
        if ( maxsize - i + 2 * k == j ) {
          val += (Real)comb(i, k) * a.col(asize - 1 - i) * pow(freqSq, (Real)(i-k)) / pow(bandwidth, (Real)i);
        }
      }
    }

    (*aout).col(maxasize - j) = val;
  }

  normalize((*bout), (*aout));
  return;
}

void normalize(MatrixXC& b, MatrixXC& a) {

  for (int i = 0; i < b.cols(); i++ ) {
    b.col(i).cwise() /= a.col(0);
  }

  for (int i = a.cols()-1; i >= 0; i-- ) {
    a.col(i).cwise() /= a.col(0);
  }
}

void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();
  const int maxsize = max(asize, bsize);
  
  const int rows = a.rows();
  
  (*aout).resize(rows, maxsize);
  (*bout).resize(rows, maxsize);
  
  MatrixXC val;
  
  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += comb(i, k) * comb(maxsize - i - 1, l) * b.col(bsize - 1 - i) * pow(2*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*bout).col(j) = val.real();
  }

  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += comb(i, k) * comb(maxsize - i - 1, l) * a.col(asize - 1 - i) * pow((Real)2.0*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*aout).col(j) = val.real();
  }
  
}

void bilinear(const MatrixXC& b, const MatrixXC& a, MatrixXR*  bout, MatrixXR*  aout) {
  bilinear(b, a, 1.0, bout, aout);
}
