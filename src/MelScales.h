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

#ifndef MELSCALES_H
#define MELSCALES_H

#include "Typedefs.h"
#include "Debug.h"

/**
 *
 * Mel scales computed using the Greenwood function:
 *
 * Greenwood, DD. (1990)
 * A cochlear frequency-position function for several species - 29 years later,
 * Journal of the Acoustical Society of America, vol. 87, pp. 2592-2605.
 *
 */
Real linearToMelGreenwood1990(Real linearFreq);

Real melToLinearGreenwood1990(Real melFreq);

void linearToMelMatrixGreenwood1990(const MatrixXR& linearFreq, MatrixXR* melFreq);

void melToLinearMatrixGreenwood1990(const MatrixXR& melFreq, MatrixXR* linearFreq);


/**
 *
 * Mel scales computed using the original formula proposed by:
 *
 * Stevens, Stanley Smith; Volkman; John; & Newman, Edwin. (1937). 
 * A scale for the measurement of the psychological magnitude of pitch.
 * Journal of the Acoustical Society of America, 8 (3), 185-190.
 *
 */
Real linearToMelStevens1937(Real linearFreq);

Real melToLinearStevens1937(Real melFreq);

void linearToMelMatrixStevens1937(const MatrixXR& linearFreq, MatrixXR* melFreq);

void melToLinearMatrixStevens1937(const MatrixXR& melFreq, MatrixXR* linearFreq);


/**
 *
 * Mel scales computed using the formula proposed by:
 *  
 * Fant, Gunnar. (1968).
 * Analysis and synthesis of speech processes.
 * In B. Malmberg (Ed.), Manual of phonetics (pp. 173-177). Amsterdam: North-Holland.
 *
 */
Real linearToMelFant1968(Real linearFreq);

Real melToLinearFant1968(Real melFreq);

void linearToMelMatrixFant1968(const MatrixXR& linearFreq, MatrixXR* melFreq);

void melToLinearMatrixFant1968(const MatrixXR& melFreq, MatrixXR* linearFreq);

#endif  /* MELSCALES_H */
