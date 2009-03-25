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

#ifndef SPECTRALNOISESUPPRESSION_H
#define SPECTRALNOISESUPPRESSION_H

#include <Eigen/StdVector>

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

class SpectralNoiseSuppression {
protected:
  int _fftSize;
  Real _samplerate;

  Real _lowFrequency;
  Real _highFrequency;
  
  int _k0;
  int _k1;

  MatrixXR _g;

  Bands _bands;

public:
  SpectralNoiseSuppression( int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 6000.0, Real samplerate = 44100.0 );

  ~SpectralNoiseSuppression();

  void setup();
  void reset();

  void process( const MatrixXR& spectrum, MatrixXR* noise, MatrixXR* result );

  /**
     Return the low frequency of the spectral whitening.
     The default is 50.0.

     @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
  */  
  Real lowFrequency() const;  

  /**
     Specifies the low @a frequency of the spectral whitening.
     The given @a frequency must be in the range of 0 to the samplerate / 2.
     
     @sa lowFrequency, highFrequency, setHighFrequency
  */
  void setLowFrequency( Real frequency, bool callSetup = true );

  /**
     Return the high frequency of the spectral whitening.
     The default is 6000.0.

     @sa lowFrequency, setLowFrequency, setHighFrequency
  */  
  Real highFrequency() const;  

  /**
     Specifies the high @a frequency of the spectral whitening.
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
     Returns the size of the FFT that has been performed for the input.
     The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const;

  /**
     Specifies the @a size of the FFT that has been performed for the input.
     The given @a size must be higher than 0.
     
     @sa fftSize()
  */
  void setFftSize( int size, bool callSetup = true );
};

#endif  /* SPECTRALNOISESUPPRESSION_H */
