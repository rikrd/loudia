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

#include "SpectralODFFlux.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralODFFlux::SpectralODFFlux(int fftSize) :
  SpectralODFBase(),
  _fftSize( fftSize )
{
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFFlux::~SpectralODFFlux() {}


void SpectralODFFlux::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFFLUX: Setting up...");

  reset();

  LOUDIA_DEBUG("SPECTRALODFFLUX: Finished set up...");
}


void SpectralODFFlux::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFFLUX: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();

  if ( rows < 2 ) {
    // Throw ValueError, it must have a minimum of 2 rows
  }

  (*odfValue).resize(rows - 1, 1);
  (*odfValue).row(0).setZero();
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Spectrum resized rows: " << rows );
    
  _spectrumAbs = fft.array().abs();
  
  (*odfValue) = (_spectrumAbs.block(1, 0, rows - 1, cols) \
                 - _spectrumAbs.block(0, 0, rows - 1, cols) \
                 ).array().clipUnder().rowwise().sum() / cols;
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Finished Processing");
}

void SpectralODFFlux::reset() {
  // Initial values
}
