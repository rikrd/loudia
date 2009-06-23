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

#ifndef FRAMECUTTER_H
#define FRAMECUTTER_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class FrameCutter
  *
  * @brief Algorithm to create frames from an input stream such as the one
  * provided by the AudioLoader.
  *
  * @author Ricard Marxer
  *
  * @sa AudioLoader
  */
class FrameCutter{
public:
protected:
  int _maxInputSize;
  int _frameSize;
  int _hopSize;
  Real _defaultValue;
  int _firstSamplePosition;

  int _indexWriter;
  int _availableToWrite;
  int _availableToRead;

  int _maxFrameCount;
  
  VectorXR _buffer;
  VectorXR _row;
  
  int read(VectorXR* frame, int release);
  int write(const VectorXR& stream);

public:
  /**
     Constructs a FrameCutter object with the specified @a inputSize, 
     @a frameSize, @a hopSize settings.
     
     @param maxInputSize maximum size of the input frames,
     must be > 0.

     @param frameSize size of the output resampled frames,
     must be > 0.

     @param hopSize the hop size of the frame cutter
  */
  FrameCutter(int maxInputSize = 1024, int frameSize = 1024, int hopSize = -1, const int firstSamplePosition = 0, Real defaultValue = 0.0);

  /**
     Destroys the algorithm and frees its resources.
  */  
  ~FrameCutter();

  /**
     Performs the frame cutting.
     
     @param stream matrix of Real values.  The number of columns of @a frames 
     must be equal to the inputSize.
     
     @param resampled pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a frames and outputSize columns. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */  
  int process(const MatrixXR& stream, MatrixXR* frames);
  
  void setup();
  void reset();

  /**
     Returns the maximum input size of the algorithm.
     
     By default it is 1024.
  */
  int maxInputSize() const { return _maxInputSize; };
  
  /**
     Specifies the maximum input @a size of the algorithm.
  */
  void setMaxInputSize( int size, bool callSetup = true ) { _maxInputSize = size; if (callSetup) setup(); };

  /**
     Returns the frame size of the algorithm.
     
     By default it is 1024.
  */
  int frameSize() const { return _frameSize; };
  
  /**
     Specifies the frame @a size of the algorithm.
  */
  void setFrameSize( int size, bool callSetup = true ) { _frameSize = size; if (callSetup) setup(); };

  /**
     Returns the hop size of the algorithm.
     
     By default it is 1024.
  */
  int hopSize() const {
    if ( _hopSize <= 0 ) {
      return std::max(_frameSize / 2, 1);
    }
    
    return _hopSize;
  };
  
  /**
     Specifies the hop @a size of the algorithm.
  */
  void setHopSize( int size, bool callSetup = true ) { _hopSize = size; if (callSetup) setup(); };

  void setFirstSamplePosition( int position, const bool callSetup = true ) { _firstSamplePosition = position; if ( callSetup ) setup(); };
  int firstSamplePosition() const { return _firstSamplePosition; };  
  
  int maxFrameCount() const;
};

#endif //#ifndef FRAMECUTTER_H
