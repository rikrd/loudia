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

#include "MelScales.h"

Real linearToMelGreenwood1990(Real linearFreq) {
  return log10((linearFreq / 165.4) + 1.0) / 2.1;
}

Real melToLinearGreenwood1990(Real melFreq) {
  return 165.4 * (pow(10.0, 2.1 * melFreq) - 1.0);
}

void linearToMelMatrixGreenwood1990(const MatrixXR& linearFreq, MatrixXR* melFreq) {
  LOUDIA_DEBUG("MELBANDS: Scaling (Greenwood 1990) linearFreq: " << linearFreq);

  (*melFreq) = ((linearFreq / 165.4).array() + 1.0).logN(10) / 2.1;
}

void melToLinearMatrixGreenwood1990(const MatrixXR& melFreq, MatrixXR* linearFreq) {
  LOUDIA_DEBUG("MELBANDS: Scaling (Greenwood 1990) melFreq: " << melFreq);

  (*linearFreq) = 165.4 * ((melFreq * 2.1).array().expN(10.0) - 1.0);
}


Real linearToMelStevens1937(Real linearFreq) {
  return log((linearFreq / 700.0) + 1.0) * 1127.01048;
}

Real melToLinearStevens1937(Real melFreq) {
  return (exp(melFreq / 1127.01048) - 1.0) * 700.0;
}

void linearToMelMatrixStevens1937(const MatrixXR& linearFreq, MatrixXR* melFreq) {
  (*melFreq) = ((linearFreq / 700.0).cwise() + 1.0).cwise().log() * 1127.01048;
}

void melToLinearMatrixStevens1937(const MatrixXR& melFreq, MatrixXR* linearFreq) {
  (*linearFreq) = ((melFreq / 1127.01048).cwise().exp().cwise() - 1.0) * 700.0;
}


Real linearToMelFant1968(Real linearFreq) {
  return (1000.0 / log(2.0)) * log(1.0 + linearFreq / 1000.0);
}

Real melToLinearFant1968(Real melFreq) {
  return 1000.0 * (exp(melFreq * log(2.0) / 1000.0) - 1.0);
}

void linearToMelMatrixFant1968(const MatrixXR& linearFreq, MatrixXR* melFreq) {
  (*melFreq) = (1000.0 / log(2.0)) * ((linearFreq / 1000.0).cwise() + 1.0).cwise().log();
}

void melToLinearMatrixFant1968(const MatrixXR& melFreq, MatrixXR* linearFreq) {
  (*linearFreq) = 1000.0 * ((melFreq * log(2.0) / 1000.0).cwise().exp().cwise() - 1.0);
}
