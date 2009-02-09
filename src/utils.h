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

#ifndef UTILS_H
#define UTILS_H

#include "typedefs.h"
#include <limits>


/**
 * Given a matrix of polynomes (one per column)
 * returns a matrix of roots (a vector of roots per column)
 */
void roots(const MatrixXR& poly, MatrixXC* result);

/**
 * Given a matrix of roots (a vector of roots per column)
 * returns a matrix of polynomes (a polynome per vector of roots)
 */
void poly(const MatrixXC& roots, MatrixXC* result);

/**
 * Given two row matrices 
 * returns the convolution of both
 */
void convolve(const MatrixXC& a, const MatrixXC& b, MatrixXC* c);
void convolve(const MatrixXR& a, const MatrixXR& b, MatrixXR* c);


/**
 * Given two row matrices 
 * returns the correlation of both
 */
void correlate(const MatrixXC& a, const MatrixXC& b, MatrixXC* c, 
               int _minlag = -std::numeric_limits<Real>::infinity(), 
               int _maxlag = std::numeric_limits<Real>::infinity());

void correlate(const MatrixXR& a, const MatrixXR& b, MatrixXR* c, 
               int _minlag = -std::numeric_limits<Real>::infinity(), 
               int _maxlag = std::numeric_limits<Real>::infinity());

/**
 * Given a row matrix 
 * returns the autocorrelation
 */
void autocorrelate(const MatrixXR& a, MatrixXR* c, 
                   int _minlag = 0, 
                   int _maxlag = std::numeric_limits<Real>::infinity());

void autocorrelate(const MatrixXC& a, MatrixXC* c, 
                   int _minlag = 0, 
                   int _maxlag = std::numeric_limits<Real>::infinity());


/**
 * Reverse in place the order of the columns
 */
void reverseCols(MatrixXC* in);
void reverseCols(MatrixXR* in);


/**
 * Calculate inplace the cumulative sum
 */
void rowCumsum(MatrixXR* in);
void colCumsum(MatrixXR* in);

/**
 * Calculate inplace shift of a matrix
 */
void rowShift(MatrixXR* in, int num);
void colShift(MatrixXR* in, int num);

/**
 * Calculate inplace range matrix
 */
void range(Real start, Real end, int steps, MatrixXC* in);
void range(Real start, Real end, int steps, int rows, MatrixXC* in);
void range(Real start, Real end, int steps, MatrixXR* in);
void range(Real start, Real end, int steps, int rows, MatrixXR* in);

/**
 * Create a matrix of complex numbers given the polar coordinates
 */
void polar(const MatrixXR& mag, const MatrixXR& phase, MatrixXC* complex);

/**
 * Convert from the b and a coefficients of an IIR filter to the
 * zeros, poles and gain of the filter
 */
void coeffsToZpk(const MatrixXR& b, const MatrixXR& a, MatrixXC* zeros, MatrixXC* poles, Real* gain);

/**
 * Convert from zeros, poles and gain of an IIR filter to the
 * the b and a coefficients of the filter
 */
void zpkToCoeffs(const MatrixXC& zeros, const MatrixXC& poles, Real gain, MatrixXC*  b, MatrixXC*  a);

/**
 * Convert from the b and a coefficients of an IIR filter with critical frequency 1.0
 * to the coefficients with the critical frequency passed as argument
 */
void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout);

/**
 * Apply the biliniear transformations to a set of coefficients
 * 
 */
void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout);

/**
 * Calculate the combinations of N elements in groups of k
 *   
 */
int comb(int N, int k);

/**
 * Calculate the aliased cardinal sine defined as:
 *   
 *   asinc(M, T, x) = sin(M * pi * x * T) / sin(pi * x * T)
 */
Real asinc(int M, Real omega);

/**
 * Calculate the Fourier transform of a hamming window
 */
void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        Real alpha, Real beta, 
                        MatrixXR* spectrum, int* begin, int* end, int bandwidth);

void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        Real alpha, Real beta, 
                        MatrixXR* spectrum, int bandwidth);

void hannTransform(Real position, Real magnitude, 
                   int windowSize, int fftSize,
                   MatrixXR* spectrum, int bandwidth = 4);

void hannTransform(Real position, Real magnitude, 
                   int windowSize, int fftSize,
                   MatrixXR* spectrum, int* begin, int* end, int bandwidth = 4);


void hammingTransform(Real position, Real magnitude, 
                      int windowSize, int fftSize,
                      MatrixXR* spectrum, int bandwidth = 4);

void hammingTransform(Real position, Real magnitude, 
                      int windowSize, int fftSize,
                      MatrixXR* spectrum, int* begin, int* end, int bandwidth = 4);

void dbToMag(const MatrixXR& db, MatrixXR* mag);

void magToDb(const MatrixXR& mag, MatrixXR* db, Real minMag = 0.0001 );

void unwrap(const MatrixXR& phases, MatrixXR* unwrapped);

void freqz(const MatrixXR& b, const MatrixXR& a, const MatrixXR& w, MatrixXC* resp);
void freqz(const MatrixXR& b, const MatrixXR& w, MatrixXC* resp);

void derivate(const MatrixXR& a, MatrixXR* b);

#endif  /* UTILS_H */
