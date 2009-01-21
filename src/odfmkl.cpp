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

#include "odfmkl.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFMKL::ODFMKL(int fftLength, Real minSpectrum) :
  ODFBase(),
  _fftLength(fftLength),
  _minSpectrum(minSpectrum)
{
  
  DEBUG("ODFMKL: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFMKL::~ODFMKL() {}


void ODFMKL::setup() {
  // Prepare the buffers
  DEBUG("ODFMKL: Setting up...");

  reset();

  DEBUG("ODFMKL: Finished set up...");
}


void ODFMKL::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFMKL: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  if ( rows < 2 ) {
    // Throw ValueError, it must have a minimum of 2 rows
  }

  (*odfValue).resize(rows - 1, 1);
  _spectrumAbs.resize(rows, halfCols);

  DEBUG("ODFMKL: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrumAbs = fft.block(0, 0, rows, halfCols).cwise().abs();

  (*odfValue) = (_spectrumAbs.block(1, 0, rows-1, cols).cwise() \
                 * (_spectrumAbs.block(1, 0, rows-1, cols).cwise() \
                    / (_spectrumAbs.block(0, 0, rows-1, cols).cwise().clipUnder(_minSpectrum))).cwise().clipUnder(_minSpectrum).cwise().logN(2.0)).rowwise().sum() / cols;
  
  DEBUG("ODFMKL: Finished Processing");
}

void ODFMKL::reset() {
  // Initial values
}
