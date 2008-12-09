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

#include "aok.h"
#include "typedefs.h"

#include <Eigen/Core>
#include <iostream>
#include <fstream>

using namespace std;

void loadFile(string filename, MatrixXR* result, int rows, int cols) {
  FILE* in = fopen( filename.c_str(), "r");
  Real coeff;
  for ( int i = 0; i<rows; i++ ) {
    for (int j = 0; j<cols; j++) {
      int r = fscanf(in, "%f", &coeff);
      (*result)(i, j) = coeff;
    }
  }
}

int main() {
  int windowSize = 64;
  int hopSize = 1;
  int fftLength = 64;
  int numFrames = 128 + 2 * windowSize - 2;
  Real normVolume = 3;
  int frameSize = 2.42 * windowSize + 3;

  MatrixXR in = MatrixXR::Zero(numFrames, frameSize);
  loadFile("/home/rmarxer/dev/ricaudio/src/tests/chirp.frames", &in, numFrames, frameSize);
  
  //cerr << in << endl;
  
  AOK aok(windowSize, hopSize, fftLength, normVolume);
  aok.setup();

  MatrixXR result(numFrames, fftLength);
  
  aok.process(in, &result);
  
  //cout << result << endl;

  return 0;
}

