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

#include "typedefs.h"
#include "debug.h"

#include "odfhfc.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFHFC::ODFHFC(int fftLength) :
  ODFBase(),
  _fftLength(fftLength)
{
  
  DEBUG("ODFHFC: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFHFC::~ODFHFC() {}


void ODFHFC::setup() {
  // Prepare the buffers
  DEBUG("ODFHFC: Setting up...");

  // Create the vector with the weights (weights are the frequency bin indices)
  const int halfFFTlen = (int)ceil(_fftLength / 2.0);
  range(0, halfFFTlen, halfFFTlen, &_freqBin);
  
  reset();

  DEBUG("ODFHFC: Finished set up...");
}


void ODFHFC::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFHFC: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  (*odfValue).resize(rows, 1);
  _spectrumAbs.resize(rows, halfCols);

  DEBUG("ODFHFC: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrumAbs = fft.block(0, 0, rows, halfCols).cwise().abs();  

  for (int row = 0; row < rows; row ++) {
    (*odfValue).row(row) = _spectrumAbs.row(row) * _freqBin.block(0, 0, 1, halfCols).transpose() / halfCols;
  }
  
  DEBUG("ODFHFC: Finished Processing");
}

void ODFHFC::reset() {
  // Initial values
}
