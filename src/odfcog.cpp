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

#include "odfcog.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFCOG::ODFCOG(int fftLength, int peakCount, int bandwidth) :
  ODFBase(),
  _fftLength(fftLength),
  _peakCount(peakCount),
  _bandwidth(bandwidth),
  _peaker(peakCount, bandwidth)
{
  
  DEBUG("ODFCOG: Constructor fftLength: " << _fftLength);
  
  setup();
}

ODFCOG::~ODFCOG() {}


void ODFCOG::setup() {
  // Prepare the buffers
  DEBUG("ODFCOG: Setting up...");

  _peaker.setup();
  
  reset();

  DEBUG("ODFCOG: Finished set up...");
}


void ODFCOG::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFCOG: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  
  (*odfValue).resize(rows, 1);

  DEBUG("ODFCOG: Processing the peaks");

  _peaker.process(fft, &_peakPos, &_peakMag, &_peakArg);

  DEBUG("ODFCOG: Spectrum resized rows: " << rows << " halfCols: " << halfCols);
  
  _spectrumAbs2 = fft.block(0, 0, rows, halfCols).cwise().abs2();
  unwrap(fft.block(0, 0, rows, halfCols).cwise().angle(), &_spectrumArg);
  derivate(_spectrumArg, &_spectrumArgDeriv);
  
  _cog = MatrixXR::Zero(_peakPos.rows(), _peakPos.cols());

  for(int row = 0; row < rows; row++) {
    for(int i = 0; i < _peakCount; i++){
      int start = max(0, (int)floor(_peakPos(row, i) - _bandwidth / 2));
      int end = min(halfCols, (int)ceil(_peakPos(row, i) + _bandwidth / 2));
      int bandwidth = end - start;

      if (_peakPos(row, i) != -1) {
        
        _cog(row, i) = ((-_spectrumArgDeriv).block(row, start, 1, bandwidth).cwise() * _spectrumAbs2.block(row, start, 1, bandwidth)).sum() / _spectrumAbs2.block(row, start, 1, bandwidth).sum();

      }
      
    }
  }

  (*odfValue) = _cog.cwise().clipUnder().rowwise().sum();
  
  DEBUG("ODFCOG: Finished Processing");
}

void ODFCOG::reset() {
  // Initial values
  _peaker.reset();
}
