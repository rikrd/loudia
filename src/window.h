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

#include "typedefs.h"
#include "debug.h"

class Window{
public:
  enum WindowType {
    RECTANGULAR = 0,
    HANN = 1,
    HANNING = 2,
    HAMMING = 3,
    COSINE = 4,
    BLACKMAN = 5,
    BLACKMANHARRIS = 6,
    NUTTALL = 7,
    BLACKMANNUTTALL = 8,
    CUSTOM = 9
  };

protected:
  int _frameSize;
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
  Window(int frameSize, WindowType windowType = RECTANGULAR);
  ~Window();
  
  void process(const MatrixXC& frames, MatrixXC* windowedFrames);
  void process(const MatrixXR& frames, MatrixXR* windowedFrames);
  void process(const MatrixXR& frames, MatrixXC* windowedFrames);
  
  void setup();
  void reset();

  int frameSize() const;
  WindowType windowType() const;
  MatrixXR window() const;
  void setWindow(MatrixXR window);
};

#endif  /* WINDOW_H */
