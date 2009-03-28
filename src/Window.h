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

#ifndef WINDOW_H
#define WINDOW_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class Window
  *
  * @brief Algorithm to create and apply several type of windows on vectors of Real 
  * or Complex values.
  *
  * This class represents an object to apply a window on frames of Real or Complex data.
  * The algorithm takes as input N-point vectors of Real (or Complex) values 
  * and returns N-point vectors of Real (or Complex) values where the samples are weighted
  * by a weighting window.
  *
  * 5 types of windows are implemented:
  * -# Rectangular
  * -# Hann or Hanning
  * -# Hamming
  * -# Cosine
  * -# Blackmann
  * -# Blackmann Harris
  * -# Nuttall
  * -# Blackman Nuttall
  *
  * The Window type can be selected using the 
  * setWindowType() method.
  *
  * Additionally a Custom window can be specified using 
  * the setWindow() method.
  *
  * @author Ricard Marxer
  *
  * @sa FFT
  */
class Window{
public:
  /**
    @enum WindowType
    @brief Specifies the type of the window.

    @sa windowType
  */
  enum WindowType {
    RECTANGULAR         = 0 /**< Rectangular window */,
    HANN                = 1 /**< Hann window */,
    HANNING             = 2 /**< Alias for a Hann window */,
    HAMMING             = 3 /**< Hamming window */,
    COSINE              = 4 /**< Cosine window */,
    BLACKMAN            = 5 /**< Blackman window */,
    BLACKMANHARRIS      = 6 /**< Blackman-Harris window */,
    NUTTALL             = 7 /**< Nuttall window */,
    BLACKMANNUTTALL     = 8 /**< Blackman-Nuttall window */,
    CUSTOM              = 9 /**< Custom window. Note that this window type must be select
                             when setting the window using setWindow()*/
  };

protected:
  int _inputSize;
  WindowType _windowType;
  MatrixXR _window;
  
  MatrixXR hann(int length);
  MatrixXR hamming(int length);
  MatrixXR cosine(int length);

  MatrixXR blackmanType(int length, Real a0, Real a1, Real a2, Real a3);
  MatrixXR blackman(int length);
  MatrixXR nuttall(int length);
  MatrixXR blackmanHarris(int length);
  MatrixXR blackmanNuttall(int length);

  template<typename FrameMatrixType, typename WindowedMatrixType>
  void process(const FrameMatrixType& frames, WindowedMatrixType* windowedFrames);
 
public: 
  /**
     Constructs a Window object with the given @a inputSize and @a windowType parameters
     given.
  */
  Window(int inputSize = 1024, WindowType windowType = RECTANGULAR);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~Window();

  void setup();
  void reset();

  /**
     Applies the window on each of the rows of @a frames and
     puts the result in the rows of @a windowedFrames.
     
     @param frames matrix of Real (or Complex) values.  The number of columns of @a frames must
     be equal to the inputSize property.
     
     @param windowedFrames pointer to a matrix of Real (or Complex) values for the output.
     The matrix should have the same number of rows and columns as @a frames. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */  
  void process(const MatrixXC& frames, MatrixXC* windowedFrames);
  void process(const MatrixXR& frames, MatrixXR* windowedFrames);
  void process(const MatrixXR& frames, MatrixXC* windowedFrames);

  /**
     Returns the input size of the algorithm.
     
     By default it is 1024.
  */
  int inputSize() const;
  
  /**
     Specifies the input size of the algorithm.
  */
  void setInputSize( int size, bool callSetup = true );
  
  /**
     Return the type of the window
     
     By default it is RECTANGULAR.
  */
  WindowType windowType() const;
  
  /**
     Specify the type of the window.
  */
  void setWindowType( WindowType type, bool callSetup = true );

  /**
     Return the single row matrix of Real values representing the window.
     
     The number of cols of the window will be equal to inputSize.
     
     By default it is a single row matrix with all values set to 1.0.
  */  
  const MatrixXR& window() const;

  /**
     Specify the single row matrix of Real values representing the window.
     
     The number of cols of the window must be equal to inputSize.
     
     Note that when the window is set, using setWindow(),
     the window type is automatically set to CUSTOM.

     By default it is a single row matrix with all values set to 1.0.
  */
  void setWindow( const MatrixXR& window, bool callSetup = true );
};

#endif  /* WINDOW_H */
