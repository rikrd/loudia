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

#ifndef FILTERUTILS_H
#define FILTERUTILS_H

#include "Typedefs.h"
#include "Debug.h"

#include "Utils.h"

/**
 * Create the zeros, poles and gain of an analog prototype of a Chebyshev Type I filter.
 */
void chebyshev1(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

/**
 * Create the  zeros, poles and gain of an analog prototype of a Chebyshev Type II filter.
 */
void chebyshev2(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

/**
 * Create the zeros, poles and gain of an analog prototype of a Butterworth filter.
 */
void butterworth(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

/**
 * Create the zeros, poles and gain of an analog prototype of a Bessel filter.
 */
void bessel(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

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
 * Convert from the b and a coefficients from low pass to low pass of an IIR filter with critical frequency 1.0
 * to the coefficients with the critical frequency passed as argument
 */
void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout);

/**
 * Convert from the b and a coefficients from low pass to high pass of an IIR filter with critical frequency 1.0
 * to the coefficients with the critical frequency passed as argument
 */
void lowPassToHighPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout);

/**
 * Convert from the b and a coefficients from low pass to band pass of an IIR filter with critical frequency 1.0
 * to the coefficients with the critical frequency passed as argument
 */
void lowPassToBandPass(const MatrixXC& b, const MatrixXC& a, Real freq, Real freqStop, MatrixXC*  bout, MatrixXC*  aout);

/**
 * Convert from the b and a coefficients from low pass to band stop of an IIR filter with critical frequency 1.0
 * to the coefficients with the critical frequency passed as argument
 */
void lowPassToBandStop(const MatrixXC& b, const MatrixXC& a, Real freq, Real freqStop, MatrixXC*  bout, MatrixXC*  aout);

/**
 * Normalize to a first coefficient
 * 
 */
void normalize(MatrixXC& b, MatrixXC& a);

/**
 * Apply the biliniear transformations to a set of coefficients
 * 
 */
void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout);

#endif  /* FILTERUTILS_H */
