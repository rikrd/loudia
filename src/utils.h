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

/**
 * Create a matrix of complex numbers given the polar coordinates
 */
void polar(MatrixXR mag, MatrixXR phase, MatrixXC* complex);

/**
 * Convert from the b and a coefficients of an IIR filter to the
 * zeros, poles and gain of the filter
 */
void coeffsToZpk(MatrixXR b, MatrixXR a, MatrixXC* zeros, MatrixXC* poles, Real* gain);

#endif  /* UTILS_H */
