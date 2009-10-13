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

#ifndef BARKBANDS_H
#define BARKBANDS_H

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

/**
  * @class BarkBands
  *
  * @brief Algorithm to get the band values of Bark-scale frequency warpped magnitude spectrums.
  *
  * This class represents an object to calculate the Bark bands on vectors
  * of Real values.  This method is a special case of the Bands algorithm.
  *
  * The method consists in a set rectangular windows spaced evenly on a
  * Bark-frequency scale.
  *
  * The sampleRate and FFT size of the input spectrum are specified using setSampleRate() and
  * setFftSize().
  *
  * The first and last bands are specified using setLowBand() and
  * setHighBand().
  *
  *
  * @author Ricard Marxer
  *
  * @sa Bands, MelBands
  */
class BarkBands {
protected:
  Real _lowBand;
  Real _highBand;
  Real _sampleRate;
  int _fftSize;

  Bands _bands;
  MatrixXR _centersLinear;

public:
  /**
     Constructs a Bark bands object with the specified @a lowBand, @a highBand,
     @a bandCount, @a sampleRate, @a fftSize and @a scaleType settings.

     @param lowBand band of the lowest Bark band,
     must be greater than zero 0 and lower than half the sampleRate.

     @param highBand band of the highest Bark band,
     must be greater than zero 0 and lower than half the sampleRate.

     @param bandCount number of Bark bands.

     @param sampleRate sampleRate frequency of the input signal.

     @param fftSize size of the FFT.
  */
  BarkBands(int lowBand = 0, int highBand = 23, Real sampleRate = 44100.0, int fftSize = 1024);

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
     Return the low band.
     The default is 0.

     @sa lowBand, highBand, setLowBand, setHighBand
  */
  Real lowBand() const;

  /**
     Specifies the low @a band of the spectral whitening.
     The given @a band must be in the range of 0 to the sampleRate / 2.

     @sa lowBand, highBand, setHighBand
  */
  void setLowBand( Real band, bool callSetup = true );

  /**
     Return the high band.
     The default is 23.

     @sa lowBand, setLowBand, setHighBand
  */
  Real highBand() const;

  /**
     Specifies the high @a band.
     The given @a band must be in the range of 0 to the sampleRate / 2.

     @sa lowBand, highBand, setLowBand
  */
  void setHighBand( Real band, bool callSetup = true );

  /**
     Return the sampleRate frequency of the input signal.
     The default is 44100.0.

     @sa setSampleRate
  */
  Real sampleRate() const;

  /**
     Specifies the sampleRate @a frequency of the input signal.

     @sa sampleRate
  */
  void setSampleRate( Real frequency, bool callSetup = true );

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

#endif  /* BARKBANDS_H */
