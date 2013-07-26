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

#include "SpectralODFHFC.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralODFHFC::SpectralODFHFC(int fftSize) :
  SpectralODFBase()
{
  
  LOUDIA_DEBUG("SPECTRALODFHFC: Constructor fftSize: " << _fftSize);
  setFftSize( fftSize, false );
  setup();
}

SpectralODFHFC::~SpectralODFHFC() {}


void SpectralODFHFC::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFHFC: Setting up...");

  SpectralODFBase::setup();

  // Create the vector with the weights (weights are the frequency bin indices)
  range(0, _halfSize, _halfSize, &_freqBin);
  
  reset();

  LOUDIA_DEBUG("SPECTRALODFHFC: Finished set up...");
}


void SpectralODFHFC::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFHFC: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();

  (*odfValue).resize(rows, 1);
  
  _spectrumAbs = fft.array().abs();

  for (int row = 0; row < rows; row ++) {
    (*odfValue).row(row) = _spectrumAbs.row(row) * _freqBin.block(0, 0, 1, cols).transpose() / cols;
  }
  
  LOUDIA_DEBUG("SPECTRALODFHFC: Finished Processing");
}

void SpectralODFHFC::reset() {
  // Initial values
}
