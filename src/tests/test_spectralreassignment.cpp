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

#include "SpectralReassignment.h"
#include "Window.h"

#include "Typedefs.h"
#include "Debug.h"

using namespace std;

int main() {
  int frameLength = 128;
  int fftLength = 256;
  Real samplerate = 8000;
  Window::WindowType windowType(Window::HAMMING);

  MatrixXR in = MatrixXR::Constant(1, frameLength, 2.0);
  MatrixXR in2(1, frameLength);
  for (int i=0; i< in2.cols(); i++){
    in2(0, i) = cos(2.0 * M_PI * 440.0 * Real(i) / samplerate);
  }
  MatrixXR in3(1, frameLength);
  for (int i=0; i< in2.cols(); i++){
    in3(0, i) = cos(2.0 * M_PI * 440.0 * Real(i) / samplerate) + cos(2.0 * M_PI * 220.0 * Real(i) / samplerate);
  }
  
  SpectralReassignment specreassign(frameLength, fftLength, samplerate, windowType);
  specreassign.setup();

  MatrixXC result(1, fftLength);
  
  specreassign.process(in, &result);
  cout << "result in:" << result.cwise().abs() << endl;

  specreassign.process(in2, &result);
  cout << "result in2:" << result.cwise().abs() << endl;

  specreassign.process(in3, &result);
  cout << "result in3:" << result.cwise().abs() << endl;

  return 0;
}

