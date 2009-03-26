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

#ifndef SPECTRALWHITENING_H
#define SPECTRALWHITENING_H

#include "Typedefs.h"
#include "Debug.h"

#include "MelBands.h"

/**
  * @class SpectralWhitening
  *
  * @brief Algorithm to whiten the magnitude of spectrums represented as vectors of Real values.
  *
  * This class represents an object to perform spectral whitening on vectors 
  * of Real values.  Which is a useful technique to make the peaks of a spectrum magnitude 
  * stand out in harmonic sounds.
  *
  * This implementation consists in calculating the Mel bands and create a linear interpolation
  * between these to be used as a wheighting parameter of the spectrum's compression. 
  *
  * The samplerate and FFT size of the input spectrum are specified using setSamplerate() and
  * setFftSize().
  *
  * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
  * setHighFrequency().
  *
  * The number of Mel bands is specified using setBandCount().
  *
  * The compression factor of the whitening process is specified by setCompressionFactor().
  *
  * @author Ricard Marxer
  *
  * @sa MelBands, Bands, PeakDetection
  */
class SpectralWhitening {
protected:
  int _fftSize;
  int _halfSize;
  Real _lowFrequency;
  Real _highFrequency;

  Real _samplerate;
  Real _compressionFactor;
  int _bandCount;

  MelBands::ScaleType _scaleType;

  MatrixXR _centers;

  MatrixXR _bandEnergy;
  MatrixXR _compressionWeights;

  MelBands _bands;

public:
  /**
     Constructs a spectral whitening object with the specified @a lowFrequency, @a highFrequency, 
     @a bandCount, @a samplerate, @a fftSize, @a compressionFactor and @a scaleType settings.
     
     @param lowFrequency frequency of the lowest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
     
     @param highFrequency frequency of the highest Mel band,
     must be greater than zero 0 and lower than half the samplerate.
 
     @param bandCount number of Mel bands.
     
     @param samplerate samplerate frequency of the input signal.

     @param fftSize size of the FFT.
     
     @param compressionFactor factor of the compression process in the whitening.
     
     @param scaleType scale used for the frequency warping.
  */  
  SpectralWhitening(int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 6000.0, Real samplerate = 44100.0, Real compressionFactor = 0.33, int bandCount = 40, MelBands::ScaleType scaleType = MelBands::GREENWOOD);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~SpectralWhitening();

  void setup();
  void reset();

  /**
     Performs a whitening on each of the rows of @a spectrums.
     Puts the resulting whitened spectrum in the rows of @a whitened.
     
     @param spectrums matrix of Real values representing one spectrum magnitude per row.
     The number of columns of @a spectrum must be equal to the fftSize / 2 + 1 where 
     fftSize is specified using setFftSize().
     
     @param whitened pointer to a matrix of Real values representing one whitened spectrum per row.
     The matrix should have the same number of rows and columns as @a spectrums.
     
     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& spectrums, MatrixXR* whitened);

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

  /**
     Return the compression factor of the whitening.
     The default is 0.33.

     @sa setCompressionFactor
  */  
  Real compressionFactor() const;  

  /**
     Specifies the compression @a factor of the whitening.
     
     @sa compressionFactor
  */
  void setCompressionFactor( Real factor, bool callSetup = true );  
  
};

#endif  /* SPECTRALWHITENING_H */
