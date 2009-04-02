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

#include "PeakCOG.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

PeakCOG::PeakCOG(int fftLength, int bandwidth) :
  _fftLength(fftLength),
  _bandwidth(bandwidth)
{
  
  DEBUG("PEAKCOG: Constructor fftLength: " << _fftLength);
  
  setup();
}

PeakCOG::~PeakCOG() {}


void PeakCOG::setup() {
  // Prepare the buffers
  DEBUG("PEAKCOG: Setting up...");
 
  reset();

  DEBUG("PEAKCOG: Finished set up...");
}


void PeakCOG::process(const MatrixXC& fft, const MatrixXR& peakPos, MatrixXR* peakCog) {
  DEBUG("PEAKCOG: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  const int peakCount = peakPos.cols();
  
  DEBUG("PEAKCOG: fft.shape " << fft.rows() << "," << fft.cols());
  _spectrumAbs2 = fft.block(0, 0, rows, halfCols).cwise().abs2();
  DEBUG("PEAKCOG: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
    
  unwrap(fft.block(0, 0, rows, halfCols).cwise().angle(), &_spectrumArg);
  derivate(_spectrumArg, &_spectrumArgDeriv);
  
  (*peakCog).resize(rows, peakCount);
  (*peakCog).setZero();

  for(int row = 0; row < rows; row++) {
    for(int i = 0; i < peakCount; i++){     
      if (peakPos(row, i) != -1) {
        int start = max(0, (int)floor(peakPos(row, i) - _bandwidth / 2));
        int end = min(halfCols, (int)ceil(peakPos(row, i) + _bandwidth / 2));

        if ( (end - start) > 0) {
          (*peakCog)(row, i) = ((-_spectrumArgDeriv).block(row, start, 1, end-start).cwise() * _spectrumAbs2.block(row, start, 1, end-start)).sum() / _spectrumAbs2.block(row, start, 1, end-start).sum();
        }
      }
      
    }
  }
  
  DEBUG("PEAKCOG: Finished Processing");
}

void PeakCOG::reset() {
  // Initial values
}
