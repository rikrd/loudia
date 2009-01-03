/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#ifndef MFCC_H
#define MFCC_H

#include "typedefs.h"
#include "debug.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "melbands.h"
#include "dct.h"



//using namespace std;

// import most common Eigen types 
//using namespace Eigen;

class MFCC {
protected:
  // Internal parameters
  Real _lowFreq;
  Real _highFreq;
  int _numBands;
  Real _samplerate;
  int _spectrumLength;

  int _numCoeffs;
  
  Real _minSpectrum;
  Real _power;
  
  // Internal variables
  MelBands _melbands;
  DCT _dct;

  MatrixXR _bands;
  MatrixXR _coeffs;

public:
  MFCC(Real lowFreq, Real highFreq, int numBands, Real samplerate, int spectrumLength, int numCoeffs, Real minSpectrum = 1e-10, Real power = 1.0);

  ~MFCC();

  void setup();

  void process(MatrixXR spectrum, MatrixXR* mfccCoeffs);

  void reset();

  int numCoeffs() const;

  Real lowFreq() const;

  Real highFreq() const;

};

#endif  /* MFCC_H */
