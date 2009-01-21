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

#include "odfspectralflux.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFSpectralFlux::ODFSpectralFlux(int fftLength) :
  ODFBase(),
  _fftLength(fftLength)
{
  
  DEBUG("ODFSpectralFlux: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFSpectralFlux::~ODFSpectralFlux() {}


void ODFSpectralFlux::setup() {
  // Prepare the buffers
  DEBUG("ODFSpectralFlux: Setting up...");

  reset();

  DEBUG("ODFSpectralFlux: Finished set up...");
}


void ODFSpectralFlux::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFSpectralFlux: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  if ( rows < 2 ) {
    // Throw ValueError, it must have a minimum of 2 rows
  }

  (*odfValue).resize(rows - 1, 1);
  (*odfValue).row(0).setZero();
  
  DEBUG("ODFSpectralFlux: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrum.set(fft.block(0, 0, rows, halfCols));
  
  _spectrumAbs.set(_spectrum.cwise().abs());
  
  (*odfValue) = (_spectrumAbs.block(1, 0, rows-1, halfCols) \
                 - _spectrumAbs.block(0, 0, rows-1, halfCols)           \
                 ).cwise().clipUnder().rowwise().sum() / halfCols;
  
  DEBUG("ODFSpectralFlux: Finished Processing");
}

void ODFSpectralFlux::reset() {
  // Initial values
}
