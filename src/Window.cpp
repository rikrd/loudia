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

#include "Typedefs.h"
#include "Debug.h"

#include "Window.h"

using namespace std;
using namespace Eigen;

Window::Window(int inputSize, Window::WindowType windowType)
{
  LOUDIA_DEBUG("WINDOW: Constructor inputSize: " << inputSize << 
        ", windowType: " << windowType);

  setInputSize( inputSize, false );
  setWindow( MatrixXR::Ones(1, inputSize), false );
  setWindowType( windowType, false );

  setup();

  LOUDIA_DEBUG("WINDOW: Constructed");
}

Window::~Window(){}

void Window::setup(){
  LOUDIA_DEBUG("WINDOW: Setting up...");
  
  switch(_windowType){

  case RECTANGULAR:
    _window = MatrixXR::Ones(1, _inputSize);
    break;

  case HANN:
  case HANNING:
    _window = hann(_inputSize);
    break;

  case HAMMING:
    _window = hamming(_inputSize);
    break;

  case COSINE:
    _window = hamming(_inputSize);
    break;
    
  case BLACKMAN:
    _window = blackman(_inputSize);
    break;

  case NUTTALL:
    _window = nuttall(_inputSize);
    break;

  case BLACKMANHARRIS:
    _window = blackmanHarris(_inputSize);
    break;

  case BLACKMANNUTTALL:
    _window = blackmanNuttall(_inputSize);
    break;
  case CUSTOM:
    break;
    
  default:
    // Throw ValueError unknown window type
    break;
  }
  
  LOUDIA_DEBUG("WINDOW: Finished set up...");
}

MatrixXR Window::hann(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = 0.5 * (1 - cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0)));
  }

  return result;
}

MatrixXR Window::hamming(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = 0.53836 - 0.46164 * cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0));
  }

  return result;
}

MatrixXR Window::cosine(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = sin((M_PI * (Real)i) / ((Real)length - 1.0));
  }

  return result;
}

MatrixXR Window::blackmanType(int length, Real a0, Real a1, Real a2, Real a3){
  MatrixXR result(1, length);

  Real pi_length_1 = M_PI / ((Real)length - 1.0);

  for(int i = 0; i < length; i++ ){
    result(0, i) = a0 \
                   - a1 * cos(2.0 * (Real)i * pi_length_1) \
                   + a2 * cos(4.0 * (Real)i * pi_length_1) \
                   + a3 * cos(6.0 * (Real)i * pi_length_1);
  }

  return result;
}

MatrixXR Window::blackman(int length){
  Real a0 = (1 - 0.16) / 2.0;
  Real a1 = 0.5;
  Real a2 = 0.16 / 2.0;
  Real a3 = 0.0;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::nuttall(int length){
  Real a0 = 0.355768;
  Real a1 = 0.487396;
  Real a2 = 0.144232;
  Real a3 = 0.012604;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::blackmanHarris(int length){
  Real a0 = 0.35875;
  Real a1 = 0.48829;
  Real a2 = 0.14128;
  Real a3 = 0.01168;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::blackmanNuttall(int length){
  Real a0 = 0.3635819;
  Real a1 = 0.4891775;
  Real a2 = 0.1365995;
  Real a3 = 0.0106411;
  
  return blackmanType(length, a0, a1, a2, a3);
}

template<typename FrameMatrixType, typename WindowedMatrixType>
void Window::process(const FrameMatrixType& frames, WindowedMatrixType* windowedFrames){
  (*windowedFrames).resize(frames.rows(), _inputSize);

  for (int i = 0; i < frames.rows(); i++){
    // Process and set
    (*windowedFrames).row(i) = (frames.row(i).cwise() * _window).template cast<typename WindowedMatrixType::Scalar>();
  }

}

void Window::process(const MatrixXC& frames, MatrixXC* windowedFrames){
  process<MatrixXC, MatrixXC>(frames, windowedFrames);
}

void Window::process(const MatrixXR& frames, MatrixXC* windowedFrames){
  process<MatrixXR, MatrixXC>(frames, windowedFrames);
}

void Window::process(const MatrixXR& frames, MatrixXR* windowedFrames){
  process<MatrixXR, MatrixXR>(frames, windowedFrames);
}


void Window::reset(){
}

int Window::inputSize() const{
  return _inputSize;
}

void Window::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}


Window::WindowType Window::windowType() const{
  return _windowType;
}

void Window::setWindowType( WindowType type, bool callSetup ) {
  _windowType = type;
  if ( callSetup ) setup();
}


const MatrixXR& Window::window() const{
  return _window;
}

void Window::setWindow( const MatrixXR& window, bool callSetup ){
  if (window.cols() != _inputSize || window.rows() != 1) {
    // Throw exception wrong window size
  }

  setWindowType(CUSTOM, false);
  _window = window;

  if ( callSetup ) setup();
}
