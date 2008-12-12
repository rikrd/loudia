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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class Window{
public:
  enum WindowType {
    RECTANGULAR,
    HANN,
    HANNING,
    HAMMING,
    COSINE,
    BLACKMAN,
    BLACKMANHARRIS,
    NUTTALL,
    BLACKMANNUTTALL
  };

protected:
  int _frameSize;
  WindowType _windowType;
  MatrixXR _window;
  
  template<class F, class W>
  void process(F frames, W* windowedFrames);

  MatrixXR hann(int length);
  MatrixXR hamming(int length);
  MatrixXR cosine(int length);

  MatrixXR blackmanType(int length, Real a0, Real a1, Real a2, Real a3);
  MatrixXR blackman(int length);
  MatrixXR nuttall(int length);
  MatrixXR blackmanHarris(int length);
  MatrixXR blackmanNuttall(int length);
 
public: 
  Window(int frameSize, WindowType windowType = RECTANGULAR);
  ~Window();
  
  void process(MatrixXC frames, MatrixXC* windowedFrames);
  void process(MatrixXR frames, MatrixXR* windowedFrames);
  void process(MatrixXR frames, MatrixXC* windowedFrames);
  
  void setup();
  void reset();

  int frameSize() const;
  WindowType windowType() const;
  MatrixXR window() const;
};

#endif  /* WINDOW_H */
