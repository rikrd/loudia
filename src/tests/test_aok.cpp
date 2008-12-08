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

void loadFile(char* filename, MatrixXR* result, int rows = 2, int cols = 64) {
  FILE* in = fopen( filename, "r");
  
  for ( int i = 0; i<rows; i++ ) {
    for (int j = 0; j<cols; j++) {
      int r = fscanf(in, "%f", &((*result)(i, j)));
    }
  }
}

int main() {
  int windowSize = 64;
  int hopSize = 1;
  int fftLength = 64;
  int numFrames = 223;
  Real normVolume = 4;

  MatrixXR in = MatrixXR::Random(numFrames, windowSize);
  loadFile("/home/rmarxer/dev/ricaudio/src/tests/papertest.frames", &in, numFrames, windowSize);
  
  cout << in;

  AOK aok(windowSize, hopSize, fftLength, normVolume);
  aok.setup();

  MatrixXR result(numFrames, fftLength);
  
  aok.process(in, &result);
  cout << result << endl;

  return 0;
}

