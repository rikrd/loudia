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

#ifndef MFCC_H
#define MFCC_H

#include "Typedefs.h"
#include "Debug.h"

#include "MelBands.h"
#include "DCT.h"

class MFCC {
protected:
  // Internal parameters
  Real _lowFrequency;
  Real _highFrequency;
  int _bandCount;
  Real _samplerate;
  int _fftSize;

  int _coefficientCount;
  
  Real _minSpectrum;
  Real _power;
  
  // Internal variables
  MelBands _melbands;
  DCT _dct;

  MatrixXR _bands;
  MatrixXR _coeffs;

public:
  MFCC(Real lowFrequency = 300.0, Real highFrequency = 16000.0, int bandCount = 40.0, Real samplerate = 44100.0, int fftSize = 1024, int coefficientCount = 13, Real minSpectrum = 1e-10, Real power = 1.0);

  ~MFCC();

  void reset();
  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* mfccCoefficients);

  /**
     Returns the number of coefficients to be calculated.
     The default is 13.
     
     @sa setCoefficientCount()
  */
  int coefficientCount() const;

  /**
     Specifies the @a count of coefficients to be calculated.
     The given @a count must be in the range between 0 and (input size - 1).
          
     @sa coefficientCount()
  */
  void setCoefficientCount( int count, bool callSetup = true );

  /**
     Returns the number of bands to be performed.
     The default is 40.
     
     @sa setBandCount()
  */
  int bandCount() const;

  /**
     Specifies the @a count of bands to be performed.
          
     @sa bandCount()
  */
  void setBandCount( int count, bool callSetup = true );

  /**
     Return the low frequency of the MFCC.
     The default is 300.0.

     @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
  */  
  Real lowFrequency() const;  

  /**
     Specifies the low @a frequency of the MFCC.
     The given @a frequency must be in the range of 0 to the samplerate / 2.
     
     @sa lowFrequency, highFrequency, setHighFrequency
  */
  void setLowFrequency( Real frequency, bool callSetup = true );

  /**
     Return the high frequency of the MFCC.
     The default is 16000.0.

     @sa lowFrequency, setLowFrequency, setHighFrequency
  */  
  Real highFrequency() const;  

  /**
     Specifies the high @a frequency of the MFCC.
     The given @a frequency must be in the range of 0 to the samplerate / 2.

     @sa lowFrequency, highFrequency, setLowFrequency
  */
  void setHighFrequency( Real frequency, bool callSetup = true );

  /**
     Return the samplerate frequency of the input signal.
     The default is 44100.0.

     @sa setSamplerate
  */  
  Real samplerate() const;  

  /**
     Specifies the samplerate @a frequency of the input signal.
     
     @sa samplerate
  */
  void setSamplerate( Real frequency, bool callSetup = true );

  /**
     Returns the size of the FFT to be performed.
     The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const;

  /**
     Specifies the @a size of the FFT to be performed.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 will perform faster.
     
     @sa fftSize()
  */
  void setFftSize( int size, bool callSetup = true );

};

#endif  /* MFCC_H */
