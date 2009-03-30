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

#ifndef PITCHACF_H
#define PITCHACF_H

#include "Typedefs.h"
#include "Debug.h"

#include "PeakDetection.h"
#include "PeakInterpolation.h"
#include "Autocorrelation.h"

/**
  * @class PitchACF
  *
  * @brief Algorithm to estimate the most prominent pitch of a vector of
  * Real values reprsenting the FFT of an audio frame using the 
  * Autocorrelation function.
  *
  * This class represents an object to estimate the most prominent pitch 
  * of a vector of Real values reprsenting the FFT of an audio frame using 
  * the Autocorrelation function.
  *
  * The algorithm performs an autocorrelation on the input spectrums
  * and finds the first peak of a set of highest magnitude candidates.
  * The index of the peak and the values of the peak will determine the
  * pitch frequency and saliency.
  * 
  * The minimum peak width can be specified using setMinimumPeakWidth().
  * 
  * The number of candidates at the peak detection stage be
  * specified using setPeakCandidateCount().
  * 
  *
  * @author Ricard Marxer
  *
  * @sa PitchACF, PitchSaliency, PitchInvp
  */
class PitchACF {
protected:
  int _fftSize;
  int _halfSize;
  int _minimumPeakWidth;
  int _peakCandidateCount;

  Real _samplerate;

  PeakDetection _peak;
  PeakInterpolation _peakInterp;
  Autocorrelation _acorr;

  MatrixXR _acorred;

public:
  /**
     Constructs an autocorrelation based pitch estimation function 
     object with the specified @a fftSize, @a samplerate, @a minPeakWidth
     and @a peakCandidateCount settings.
     
     @param fftSize size of the input FFT frames, must be > 0.
     
     @param samplerate the samplerate of the input signal.  By default
     it is 1.0, so the pitches will be returned in normalized frequencies.

     @param minimumPeakWidth the minimum width of a peak in the autocorrelation
     function for it to be detected.

     @param peakCandidateCount the number of highest magnitude candidates 
     to be considered during the peak detection process of the 
     autocorrelation function.
  */
  PitchACF(int fftSize = 1024, Real samplerate = 1.0, int minimumPeakWidth = 6, int peakCandidateCount = 10);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~PitchACF();

  void setup();
  void reset();

  /**
     Performs a pitch estimation on each of the rows of @a spectrums.
     Puts the resulting estimated pitches in the rows of @a pitches
     and the saliencies of each pitch in the rows of @a saliencies.
     
     @param spectrums matrix of Real values representing one 
     spectrum magnitude per row.
     The number of columns of @a spectrum must be equal 
     to the fftSize / 2 + 1 where 
     fftSize is specified using setFftSize().
     
     @param pitches pointer to a matrix of Real values representing 
     the frequencies of the estimated pitches as rows.
     The matrix should have the same number of rows as @a spectrums and as many
     columns as the count of estimated pitches.  Note that this algorithm is
     only capable of detecting a single pitch at each frame, therefore @a pitches
     will be a single column matrix.

     @param saliencies pointer to a matrix of Real values representing
     the saliencies of the estimated pitches as rows.
     The matrix should have the same number of rows as @a spectrums and as many
     columns as the count of estimated pitches.  Note that this algorithm is
     only capable of detecting a single pitch at each frame, therefore @a pitches
     will be a single column matrix.
     
     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& spectrums, MatrixXR* pitches, MatrixXR* saliencies);
  
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
     Returns the minimum width for a peak to be detected in the
     autocorrelation function.
     The default is 6.
     
     @sa setMinimumPeakWidth()
  */
  int minimumPeakWidth() const;
  
  /**
     Specifies the minimum @a width for a peak to be detected in the
     autocorrelation function.
     
     @sa minimumPeakWidth()
  */
  void setMinimumPeakWidth( int width, bool callSetup = true );

  /**
     Returns the number of highest magnitude candidates to be considered 
     during the peak detection process of the autocorrelation function.

     Note that if the value is <= 0, then no preselection is performed
     and all detected peaks are considered as candidates.

     By default it is 6.
  */
  int peakCandidateCount() const;

  /**
     Specifies the number of highest magnitude candidates to be considered 
     during the peak detection process of the autocorrelation function.

     Note that if the value is <= 0, then no preselection is performed
     and all detected peaks are considered as candidates.
  */
  void setPeakCandidateCount( int count, bool callSetup = true );

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

};

#endif  /* PITCHACF_H */
