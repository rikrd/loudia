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

#ifndef PITCHINVERSEPROBLEM_H
#define PITCHINVERSEPROBLEM_H

#include "Typedefs.h"
#include "Debug.h"

#include "PeakDetection.h"
#include "PeakInterpolation.h"

#include <Eigen/LU>

class PitchInverseProblem {
protected:
  int _fftSize;
  int _halfSize;
  Real _lowFrequency;
  Real _highFrequency;
  int _pitchCount;
  int _harmonicCount;
  int _frequencyCandidateCount;
  Real _peakWidth;

  Real _samplerate;

  Real _tMin;
  Real _tMax;
  Real _alpha;
  Real _beta;
  Real _inharmonicity;
  Real _regularisation;

  MatrixXR _projectionMatrix;
  MatrixXR _inverseProjectionMatrix;

  PeakDetection _peak;
  PeakInterpolation _peakInterp;

  void harmonicWeight(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
  void harmonicSpread(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
  void harmonicPosition(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
  Real harmonicWeight(Real f, Real fMin, Real fMax, int harmonicIndex);
  Real harmonicSpread(Real f, Real fMin, Real fMax, int harmonicIndex);
  Real harmonicPosition(Real f, Real fMin, Real fMax, int harmonicIndex);

public:
  PitchInverseProblem(int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 2100.0, Real samplerate = 44100.0, int pitchCount = 5, int harmonicCount = 10, int frequencyCandidateCount = -1, Real peakWidth = 4);

  ~PitchInverseProblem();

  void reset();
  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* frequencies);

  void projectionMatrix(MatrixXR* matrix) const;

  /**
     Return the lowest frequency candidate.
     The default is 50.0.

     @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
  */  
  Real lowFrequency() const;  

  /**
     Specifies the lowest @a frequency candidate.
     The given @a frequency must be in the range of 0 to the samplerate / 2.
     
     @sa lowFrequency, highFrequency, setHighFrequency
  */
  void setLowFrequency( Real frequency, bool callSetup = true );

  /**
     Return the highest frequency candidate.
     The default is 2100.0.

     @sa lowFrequency, setLowFrequency, setHighFrequency
  */  
  Real highFrequency() const;  

  /**
     Specifies the highest @a frequency candidate.
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
     Returns the count of candidate frequencies in with which to discretize
     the frequency space.

     Note that if the value is <= 0, then fftSize / 2 + 1 is used.

     By default it is 6.

     @sa setFrequencyCandidateCount()
  */
  int frequencyCandidateCount() const;

  /**
     Specifies the number of highest magnitude candidates to be considered 
     during the frequency detection process of the autocorrelation function.

     Note that if the value is <= 0, then no preselection is performed
     and all detected frequencys are considered as candidates.

     @sa frequencyCandidateCount()
  */
  void setFrequencyCandidateCount( int count, bool callSetup = true );

  /**
     Returns the width of the harmonic peaks.
     The default is 8.
     
     @sa setPeakWidth()
  */
  int peakWidth() const;
  
  /**
     Specifies the @a width of the harmonic peaks.
     
     @sa peakWidth()
  */
  void setPeakWidth( int width, bool callSetup = true );

  /**
     Returns the maximum count of pitches to be estimated.

     By default it is 5.

     @sa setPitchCount()
  */
  int pitchCount() const;

  /**
     Specifies the maximum @a count of pitches to be estimated.

     @sa pitchCount()
  */
  void setPitchCount( int count, bool callSetup = true );

  /**
     Returns the maximum count of harmonics to be rendered in the projection matrix.

     By default it is 10.

     @sa setHarmonicCount()
  */
  int harmonicCount() const;

  /**
     Specifies the @a count of harmonics to be rendered in the projection matrix.

     @sa harmonicCount()
  */
  void setHarmonicCount( int count, bool callSetup = true );

};

#endif  /* PITCHINVERSEPROBLEM_H */
