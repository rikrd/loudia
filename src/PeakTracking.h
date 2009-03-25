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

#ifndef PEAKTRACKING_H
#define PEAKTRACKING_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class PeakTracking
  *
  * @brief Algorithm to find peak trajectories in vectors of Complex values representing FFT, given peak
  * positions and peak magnitudes.
  *
  * The algorithm finds a maximum number of peak trajectories and returns 
  * the positions of the trajectory positions (in fractional index units) and the trajectory magnitudes (in decibel units) in
  * separate matrices.
  * 
  * The maximum number of trajectories can be specified using setTrajectoryCount().
  * 
  * The algorithm operates by matching the peaks in the current frames to the existing trajectories.
  * During the matching process a maximum frequency change of a peak can be specified using setMaximumFrequencyChange().
  *
  * The matching process also requires a trajectory to stay unmatched during a given number of frames for the trajectory to
  * disappear and leave a slot for another trajectory to be found.  The number fo silent frames can be specified using
  * silentFrameCount().
  *
  * @author Ricard Marxer
  *
  * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
  */
class PeakTracking {
protected:
  // Internal parameters
  int _trajectoryCount;
  Real _maximumFrequencyChange;
  int _silentFrameCount;

  // Internal variables
  MatrixXR _trajPositions, _trajMagnitudes;
  MatrixXR _pastTrajPositions, _pastTrajMagnitudes;
  
  bool createTrajectory(Real peakPos, Real peakMag,
                        MatrixXR* pastTrajPositions, MatrixXR* pastTrajMagnitudes,
                        MatrixXR* trajPositions, MatrixXR* trajMagnitudes,
                        int row);
    

public:
  /**
     Constructs a peak tracking object with the given @a trajectoryCount, @a maximumFrequencyChange and @a silentFrameCount settings.
  */
  PeakTracking(int trajectoryCount = 20, Real maximumFrequencyChange = 3.0, int silentFrameCount = 3);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~PeakTracking();

  void setup();
  void reset();

  /**
     Tracks peaks on each of the rows of @a ffts, @a peakPositions, @a peakMagnitudes and
     puts the resulting trajectory indices and magnitudes in the rows of @a trajectoryPositions and 
     @a trajectoryMagnitudes respectively.
     
     @param ffts matrix of Complex values representing the FFT frames.
     
     @param peakPositions matrix of Real values for the peak positions (in fractional index units).
     The matrix should have the same number of rows as @a ffts. 

     @param peakMagnitudes matrix of Real values for the peak magnitudes (in decibel units).
     The matrix should have the same number of rows as @a ffts. 

     @param trajectoryPositions pointer to a matrix of Real values for the trajectory positions (in fractional index units).
     The matrix should have the same number of rows as @a ffts and trajectoryCount columns. 

     @param trajectoryMagnitudes pointer to a matrix of Real values for the trajectory magnitudes (in decibel units).
     The matrix should have the same number of rows as @a ffts and trajectoryCount columns. 

     Note that if the count of trajectories detected is lower than trajectoryCount some values
     of the resulting arrays will be set to -1.0 in order to indicate that it is not
     a trajectory.

     Note that if the output matrices are not of the required size they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXC& ffts, 
               const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
               MatrixXR* trajectoryPositions, MatrixXR* trajectoryMagnitudes);
  
  /**
     Returns the maximum number of trajectories to be detected by the algorithm.
     
     By default it is 20.
  */
  int trajectoryCount() const;

  /**
     Specifies the maximum trajectory @a count to be detected by the algorithm.
     If <= 0, then all possible trajectories are detected.

     @param count the maximum number of trajectories to be tracked

     @param callSetup a flag specifying whether the setup() method must be call after setting the parameter.
  */
  void setTrajectoryCount( int count, bool callSetup = true );
  
  /**
     Returns the maximum frequency change of a peak for it to be matched to an existing trajectory.
     
     The change is specified in fractional index units.
     
     By default it is 3.0.
  */
  Real maximumFrequencyChange() const;

  /**
     Specifies the maximum frequency change of a peak for it to be matched to an existing trajectory.
     
     The change is specified in fractional index units.

     @param change the maximum changed allowed between a peak and an existing trajectory

     @param callSetup a flag specifying whether the setup() method must be call after setting the parameter.
  */
  void setMaximumFrequencyChange( Real change, bool callSetup = true );

  /**
     Returns the count of frames a trajectory must stay unmatched for it
     to disappear and leave the slot for another possible trajectory.
     
     By default it is 3.
  */
  int silentFrameCount() const;

  /**
     Specifies the @a count of frames a trajectory must stay unmatched for it
     to disappear and leave the slot for another possible trajectory.
     
     @param count the number of silent frames

     @param callSetup a flag specifying whether the setup() method must be call after setting the parameter.
  */
  void setSilentFrameCount( int count, bool callSetup = true );

};

#endif  /* PEAKTRACKING_H */
