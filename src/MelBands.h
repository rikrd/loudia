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

#ifndef MELBANDS_H
#define MELBANDS_H

#include <Eigen/StdVector>

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

/**
  * @class MelBands
  *
  * @brief Algorithm to get the band values of Mel-scale frequency warpped magnitude spectrums.
  *
  * This class represents an object to Mel bands on vectors 
  * of Real values.  This method is a special case of the Bands algorithm and is used
  * in many other spectral algorithms such as MFCC and SpectralWhitening.
  *
  * The method consists in a set triangular 50% overlapping windows spaced evenly on a
  * Mel-frequency scale.
  *
  * The samplerate and FFT size of the input spectrum are specified using setSamplerate() and
  * setFftSize().
  *
  * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
  * setHighFrequency().
  *
  * The number of bands is specified using setBandCount().
  *
  * @author Ricard Marxer
  *
  * @sa Bands, MFCC, SpectralWhitening
  */
class MelBands {
public:
  /**
    @enum ScaleType
    @brief Specifies the type of the scale used.

    @sa scaleType
  */
  enum ScaleType {
    STEVENS = 0 /**< Mel scales computed using the original formula proposed by:
                 *
                 * Stevens, Stanley Smith; Volkman; John; & Newman, Edwin. (1937). 
                 * A scale for the measurement of the psychological magnitude of pitch.
                 * Journal of the Acoustical Society of America, 8 (3), 185-190.
                 *
                 */,
    FANT = 1 /**< Mel scales computed using the formula proposed by:
              *  
              * Fant, Gunnar. (1968).
              * Analysis and synthesis of speech processes.
              * In B. Malmberg (Ed.), Manual of phonetics (pp. 173-177). Amsterdam: North-Holland.
              *
              */,
    GREENWOOD = 2 /**< Mel scales computed using the Greenwood function:
                   *
                   * Greenwood, DD. (1990)
                   * A cochlear frequency-position function for several species - 29 years later,
                   * Journal of the Acoustical Society of America, vol. 87, pp. 2592-2605.
                   *
                   */
  };

protected:
  Real _lowFrequency;
  Real _highFrequency;
  int _bandCount;
  Real _samplerate;
  int _fftSize;
  ScaleType _scaleType;

  Bands _bands;
  MatrixXR _centersLinear;

  Real (*_linearToMel)(Real linearFreq);
  
  Real (*_melToLinear)(Real melFreq);
  
  void (*_linearToMelMatrix)(const MatrixXR& linearFreq, MatrixXR* melFreq);
  
  void (*_melToLinearMatrix)(const MatrixXR& melFreq, MatrixXR* linearFreq);

  void triangleWindow(MatrixXR* window, Real start, Real stop, Real center = -1, Real height = Real(1.0));
  
public:
  /**
     Constructs a Mel bands object with the specified @a lowFrequency, @a highFrequency, 
     @a bandCount, @a samplerate, @a fftSize and @a scaleType settings.
     
     @param lowFrequency frequency of the lowest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
     
     @param highFrequency frequency of the highest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
 
     @param bandCount number of Mel bands.
     
     @param samplerate samplerate frequency of the input signal.

     @param fftSize size of the FFT.
          
     @param scaleType scale used for the frequency warping.
  */
  MelBands(Real lowFrequency = 50.0, Real highFrequency = 6000.0, int bandCount = 40, Real samplerate = 44100.0, int fftSize = 1024, ScaleType scaleType = GREENWOOD);

  void setup();
  void reset();

  /**
     Calculates the bands of @a spectrums.
     
     @param spectrums matrix of Real values.
     
     @param bands pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a spectrums and as many columns as the number of bands
     as specified by bandCount. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& spectrums, MatrixXR* bands);  

  /**
     Return the vector of weights.
     
     @sa starts, bandWeights, setStartsWeights
  */
  std::vector<MatrixXR> weights() const;

  /**
     Return in @a bandWeights the weights of the band given by the index @a band.
     
     @sa weights
  */
  void bandWeights(int band, MatrixXR* bandWeights) const;

  /**
     Return in @a result the single column matrix of start indices of the bands.
  */
  void starts(MatrixXI* result) const;

  /**
     Return number of bands.
  */
  int bands() const;

  /**
     Return in @a result the single column matrix of center fractional indices of the bands.
  */
  void centers(MatrixXR* result) const;

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

  /**
     Return the type of the frequency warping scale.
     
     By default it is GREENWOOD.
     
     @sa setScaleType()
  */
  MelBands::ScaleType scaleType() const;
  
  /**
     Specify the type of the frequency warping scale.
     
     @sa scaleType()
  */
  void setScaleType( MelBands::ScaleType type, bool callSetup = true );

};

#endif  /* MELBANDS_H */
