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

#ifndef SPECTRALREASSIGNMENT_H
#define SPECTRALREASSIGNMENT_H

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "window.h"
#include "fft.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class SpectralReassignment{
protected:
  int _frameSize;
  int _fftSize;
  Real _samplerate;
  Window::WindowType _windowType;
  
  Window _windowAlgo;
  Window _windowIntegAlgo;
  Window _windowDerivAlgo;

  FFT _fftAlgo;

  MatrixXC _window;
  MatrixXC _windowInteg;
  MatrixXC _windowDeriv;

  MatrixXC _fft;
  MatrixXC _fftAbs2;
  MatrixXC _fftInteg;
  MatrixXC _fftDeriv;

  MatrixXR _time;
  MatrixXR _freq;
  MatrixXR _reassignTime;
  MatrixXR _reassignFreq;

  template<class F, class W>
  void process(F frames, W* reassigned, W* fft);
 
public: 
  SpectralReassignment(int frameSize, int fftSize, Real samplerate, Window::WindowType windowType = Window::RECTANGULAR);
  ~SpectralReassignment();
  
  void process(MatrixXC frames, MatrixXC* reassigned, MatrixXC* fft);
  void process(MatrixXR frames, MatrixXC* reassigned, MatrixXC* fft);

  void process(MatrixXC frames, MatrixXC* reassigned);
  void process(MatrixXR frames, MatrixXC* reassigned);
  
  void setup();
  void reset();

  int frameSize() const;
  int fftSize() const;

  Window::WindowType windowType() const;
};

#endif  /* SPECTRALREASSIGNMENT_H */
