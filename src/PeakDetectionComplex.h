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

#ifndef PEAKDETECTIONCOMPLEX_H
#define PEAKDETECTIONCOMPLEX_H

#include "Typedefs.h"
#include "Debug.h"

#include <limits>

/**
  * @class PeakDetectionComplex
  *
  * @brief Algorithm to find peaks in a vector of Complex values.
  *
  * This class represents an object to find peaks in a vector of Complex values.
  * The algorithm finds a maximum number of peaks and returns 
  * the indices of the peaks and the values of the peaks in
  * separate matrices.
  * 
  * The maximum number of peaks can be specified using setPeakCount().
  * 
  * The resulting peak arrays may be sorted by position or by magnitude. This can be
  * specified using setSortMethod().
  *
  * When sorting by position it may be interesting to specify a number of candidates, in order
  * to perform a preselection of the highest valued peaks before sorting.  This can be specified
  * using setCandidateCount
  *
  * The implementation consists in running a sliding windows along the vector in search of 
  * indices which whose value is the maximum of the window.  The size of the window
  * defines the minimum width of the peak.
  * 
  *
  * @author Ricard Marxer
  *
  * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
  */
class PeakDetectionComplex {
public:
  /**
    @enum SortMethod
    @brief Specifies the way to sort the peak candidates before returning them.

    @sa sortMethod
  */
  enum SortMethod {
    NONE              = 0 /**< No sorting is performed */,
    BYMAGNITUDE       = 1 /**< Sorts the peak candidates by decreasing order of magnitude */,
    BYPOSITION        = 2 /**< Sorts the peak candidates by increasing order of position */
  };

protected:
  // Internal parameters
  int _peakCount;
  int _minimumPeakWidth;
  int _candidateCount;
  Real _minimumPeakContrast;
    
  SortMethod _sortMethod;

  // Internal variables
  MatrixXR _magnitudes;
  MatrixXR _phases;

public:
  /**
     Constructs a peak detection object with the given @a peakCount, @a sort method, @a minimumPeakWidth, @a candidateCount and @a minimumPeakContrast parameters given.
  */
  PeakDetectionComplex(int peakCount = 1024 / 3, SortMethod sort = BYMAGNITUDE, int minimumPeakWidth = 3, int candidateCount = -1, Real minimumPeakContrast = 0);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~PeakDetectionComplex();

  void reset();
  void setup();
  
  /**
     Detects peaks on each of the rows of @a frames and
     puts the resulting peak indices and magnitudes in the rows of @a peakPositions and 
     @a peakMagnitudes respectively.
     
     @param frames matrix of Complex values.
     
     @param peakPositions pointer to a matrix of Real values (but always Integers) for the peak indices.
     The matrix should have the same number of rows as @a frames and peakCount columns. 

     @param peakMagnitudes pointer to a matrix of Real values for the peak magnitudes.
     The matrix should have the same number of rows as @a frames and peakCount columns. 

     @param peakPhases pointer to a matrix of Real values for the peak phases.
     The matrix should have the same number of rows as @a frames and peakCount columns. 

     Note that if the count of peaks detect is lower than peakCount some values
     of the resulting arrays will be set to -1.0 in order to indicate that it is not
     a peak.

     Note that if the output matrices are not of the required size they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXC& frames,
               MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases);

  /**
     Returns the maximum number of peaks to be detected by the algorithm.
     
     By default it is 1024 / 3.
  */
  int peakCount() const;

  /**
     Specifies the maximum peak @a count to be detected by the algorithm.
     If <= 0, then all possible peaks are detected.
  */
  void setPeakCount( int count, bool callSetup = true );

  /**
     Returns the minimum width of a peak for it to be detected.
     
     By default it is 3.
  */
  int minimumPeakWidth() const;

  /**
     Specifies the minimum @a width of a peak for it to be detected.
  */
  void setMinimumPeakWidth( int width, bool callSetup = true );

  /**
     Returns the number of highest value candidates to be considered before sorting.

     Note that if the value is <= 0, then no preselection is performed
     and all detected peaks are considered as candidates.

     By default it is -1.
  */
  int candidateCount() const;

  /**
     Specifies the number of highest value candidates to be considered before sorting.

     Note that if the value is <= 0, then no preselection is performed
     and all detected peaks are considered as candidates.
  */
  void setCandidateCount( int count, bool callSetup = true );

  /**
     Returns the minimum contrast of a peak for it to be detected.
     
     The contrast is considered of a peak is the maximum value minus the minimum value
     of all the points in the peak detection running window.
     
     By default it is 0.0.
  */
  int minimumPeakContrast() const;

  /**
     Specifies the minimum contrast of a peak for it to be detected.
     
     The contrast is considered of a peak is the maximum value minus the minimum value
     of all the points in the peak detection running window.
  */
  void setMinimumPeakContrast( Real contrast, bool callSetup = true );

  /**
     Returns the method for sorting the peaks.
     
     By default it is BYMAGNITUDE.
  */
  SortMethod sortMethod() const;

  /**
     Specifies the method for sorting the peaks.
  */
  void setSortMethod( SortMethod method, bool callSetup = true );
};

#endif  /* PEAKDETECTION_H */
