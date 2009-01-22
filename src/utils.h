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


/**
 * Given a matrix of polynomes (one per column)
 * returns a matrix of roots (a vector of roots per column)
 */
void roots(const MatrixXR& poly, MatrixXC* result);


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
 * Calculate inplace range matrix
 */
void range(Real start, Real end, int steps, MatrixXR* in);

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

#endif  /* UTILS_H */
