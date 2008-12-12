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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "window.h"

#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

Window::Window(int frameSize, Window::WindowType windowType) {
  DEBUG("WINDOW: Constructor frameSize: " << frameSize << ", windowType: " << windowType);

  _frameSize = frameSize;
  _windowType = windowType;

  DEBUG("WINDOW: Constructed");
}

Window::~Window(){
}

void Window::setup(){
  DEBUG("WINDOW: Setting up...");
  
  switch(_windowType){
  case NONE:
    _window.set(MatrixXR::Ones(1, _frameSize));
    break;
  default:
    // Throw assertion unknown window type
    break;
  }
  
  DEBUG("WINDOW: Finished set up...");
}


template<class F, class W>
void Window::process(F frames, W* windowedFrames){
  for (int i = 0; i < frames.rows(); i++){
    // Process and set
    (*windowedFrames).row(i) = frames.cwise() * _window;
  }
}

void Window::process(MatrixXC frames, MatrixXC* windowedFrames){
  process<MatrixXC, MatrixXC>(frames, windowedFrames);
}

void Window::process(MatrixXR frames, MatrixXC* windowedFrames){
  process<MatrixXR, MatrixXC>(frames, windowedFrames);
}

void Window::process(MatrixXR frames, MatrixXR* windowedFrames){
  process<MatrixXR, MatrixXR>(frames, windowedFrames);
}


void Window::reset(){
}

int Window::frameSize() const{
  return _frameSize;
}

Window::WindowType Window::windowType() const{
  return _windowType;
}

MatrixXR Window::window() const{
  return _window;
}
