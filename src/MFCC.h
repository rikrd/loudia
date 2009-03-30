/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

/**
  * @class MFCC
  *
  * @brief Algorithm to calculate the Mel-frequency Cepstrum Coefficients of vectors of Real values.
  *
  * This class represents an object to perform a 
  * Mel-frequency Cepstrum Coefficients (MFCC) on vectors 
  * of Real values.  Which is a useful technique creating a sparse representation
  * of a spectrum magnitude.  The algorithm estimates a set of M coefficients
  * which should be similar for perceptually similar sounds.
  *
  * The MFCCs are calculated by taking the Discrete Cosine Transform (DCT) 
  * of the values of the Mel log powers.  The Mel log powers are calculated 
  * by applying the logarithm to the vaules of the Mel Bands (MelBands). The Mel Bands
  * are triangular overlapping windows applied on the power spectrum mapped onto the
  * Mel-frequency scale.
  *
  * The samplerate and FFT size of the input spectrum are specified using setSamplerate() and
  * setFftSize().
  *
  * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
  * setHighFrequency().
  *
  * The number of Mel bands is specified using setBandCount().
  *
  * The number of resulting DCT coefficients is specified by setCoefficientCount().
  *
  * @author Ricard Marxer
  *
  * @sa MelBands, Bands, LPC
  */
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
  /**
     Constructs an MFCC object with the specified @a lowFrequency, @a highFrequency, 
     @a bandCount, @a samplerate, @a fftSize, @a coefficientCount, @a minSpectrum and
     @a power settings.
     
     @param lowFrequency frequency of the lowest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
     
     @param highFrequency frequency of the highest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
 
     @param bandCount number of Mel bands.
     
     @param samplerate samplerate frequency of the input signal.

     @param fftSize size of the FFT.
     
     @param coefficientCount number of DCT coefficients to be estimated.
     
     @param minSpectrum value to which the spectrum is clipped before performing the logarithm.
     
     @param power value to which to power the band values before performing the DCT.  
  */  
  MFCC(Real lowFrequency = 300.0, Real highFrequency = 16000.0, int bandCount = 40.0, Real samplerate = 44100.0, int fftSize = 1024, int coefficientCount = 13, Real minSpectrum = 1e-10, Real power = 1.0);
  
  /**
     Destroys the algorithm and frees its resources.
  */
  ~MFCC();

  void reset();
  void setup();

  /**
     Performs an MFCC on each of the rows of @a frames.
     Puts the resulting MFCC coefficients in the rows of @a mfccCoefficients.
     
     @param spectrums matrix of Real values representing one spectrum magnitude per row.
     The number of columns of @a spectrum must be equal to the fftSize / 2 + 1 where 
     fftSize is specified using setFftSize().
     
     @param mfccCoefficients pointer to a matrix of Real values for the MFCC coefficients.
     The matrix should have the same number of rows as @a spectrums and coefficientCount columns.
     
     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& spectrums, MatrixXR* mfccCoefficients);

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

  /**
     Return the factor to which the bands are powered.
     The default is 1.0.

     @sa setPower
  */  
  Real power() const;  

  /**
     Specifies the @a factor to which the bands are powered.
     
     @sa power
  */
  void setPower( Real factor, bool callSetup = true );

};

#endif  /* MFCC_H */
