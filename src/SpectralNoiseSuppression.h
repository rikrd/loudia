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

#ifndef SPECTRALNOISESUPPRESSION_H
#define SPECTRALNOISESUPPRESSION_H

#include <Eigen/StdVector>

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

/**
  * @class SpectralNoiseSuppression
  *
  * @brief Algorithm to remove the non-harmonic part of the spectrums magnitudes represented as vectors of Real values.
  *
  * This class represents an object to perform spectral noise suppresion on vectors 
  * of Real values.  Which is a useful technique to keep only the peaks of a spectrum magnitude 
  * in harmonic sounds.
  *
  * This implementation consists in estimating the spectral noise by performing a 
  * moving average on the power warped spectrum magnitude using varying bandwidths.
  * The spectral noise is then removed from the original spectrum, clipping the result to zero to avoid
  * negative values. 
  *
  * The samplerate and FFT size of the input spectrum are specified using setSamplerate() and
  * setFftSize().
  *
  * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
  * setHighFrequency().
  *
  * @author Ricard Marxer
  *
  * @sa MelBands, Bands, PeakDetection
  */
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
  /**
     Constructs a spectral noise suppression object with the specified @a lowFrequency, @a highFrequency, 
     @a samplerate and @a fftSize settings.
     
     @param lowFrequency low frequency used for the magnitude warping function,
     must be greater than zero 0 and lower than half the samplerate.
     
     @param highFrequency high frequency used for the magnitude warping function,
     must be greater than zero 0 and lower than half the samplerate.
 
     @param samplerate samplerate frequency of the input signal.

     @param fftSize size of the FFT.
  
  */
  SpectralNoiseSuppression( int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 6000.0, Real samplerate = 44100.0 );

  /**
     Destroys the algorithm and frees its resources.
  */
  ~SpectralNoiseSuppression();

  void setup();
  void reset();

  /**
     Performs the estimation and suppression of the noise on each of the rows of @a spectrums.
     Puts the resulting noise spectrums and noise-suppressed spectrums in the rows of @a whitened.
     
     @param spectrums matrix of Real values representing one spectrum magnitude per row.
     The number of columns of @a spectrum must be equal to the fftSize / 2 + 1 where 
     fftSize is specified using setFftSize().
     
     @param noises pointer to a matrix of Real values representing one noise spectrum per row.
     The matrix should have the same number of rows and columns as @a spectrums.

     @param suppressed pointer to a matrix of Real values representing one noise-suppressed spectrum per row.
     The matrix should have the same number of rows and columns as @a spectrums.
     
     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process( const MatrixXR& spectrums, MatrixXR* noises, MatrixXR* suppressed );

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
