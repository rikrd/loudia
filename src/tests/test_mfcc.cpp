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

#include "mfcc.h"


#include "typedefs.h"
#include "debug.h"

#include <Eigen/Core>
#include <iostream>

using namespace std;

int main() {
  Real lowFreq = 133.33;
  Real highFreq = 4000.0;
  int nBands = 34;
  Real samplerate = 8000.0;
  int spectrumLength = 1024;
  
  int nCoeffs = 13;
  
  MatrixXR in = MatrixXR::Constant(1, spectrumLength, 1.0);
  
  MFCC mfcc(lowFreq, highFreq, nBands, samplerate, spectrumLength, nCoeffs);
  mfcc.setup();

  MatrixXR result(1, nCoeffs);
  
  for (int i=0; i<2; i++) {   
    mfcc.process(in, &result);
    cout << "result:" << result << endl;
  }

  return 0;
}

