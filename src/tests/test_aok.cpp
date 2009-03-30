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

#include "Typedefs.h"
#include "Debug.h"

#include "AOK.h"

#include <fstream>

using namespace std;

void loadFile(string filename, MatrixXC* result, int rows, int cols) {
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
  int windowSize = 256;
  int hopSize = 128;
  int fftLength = 256;
  int numFrames = 3442;
  Real normVolume = 3;
  
  //cerr << in << endl;
  
  AOK aok(windowSize, hopSize, fftLength, normVolume);
  aok.setup();

  int frameSize = aok.frameSize();
  MatrixXC in = MatrixXC::Zero(numFrames, frameSize);
  loadFile("/home/rmarxer/dev/loudia/src/tests/test.frames", &in, numFrames, frameSize);

  MatrixXR result(numFrames, fftLength);
  
  aok.process(in, &result);
  
  cout << result << endl;

  return 0;
}

