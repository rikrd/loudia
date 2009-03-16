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

#include "ODFCOG.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

ODFCOG::ODFCOG(int fftSize, int peakCount, int bandwidth) :
  ODFBase(),
  _fftSize( fftSize ),
  _peakCount( peakCount ),
  _bandwidth( bandwidth ),
  _peaker( peakCount, PeakDetect::BYMAGNITUDE, bandwidth ),
  _peakCoger( fftSize, bandwidth )
{
  
  DEBUG("ODFCOG: Constructor fftSize: " << _fftSize);
  
  setup();
}

ODFCOG::~ODFCOG() {}


void ODFCOG::setup() {
  // Prepare the buffers
  DEBUG("ODFCOG: Setting up...");

  _peaker.setup();
  _peakCoger.setup();
  
  reset();

  DEBUG("ODFCOG: Finished set up...");
}


void ODFCOG::process(const MatrixXC& fft, MatrixXR* odfValue) {
  DEBUG("ODFCOG: Processing windowed");
  const int rows = fft.rows();
  
  (*odfValue).resize(rows, 1);

  DEBUG("ODFCOG: Processing the peaks");

  _peaker.process(fft.cwise().abs(), &_peakPos, &_peakMag);

  _peakCoger.process(fft, _peakPos, &_cog);

  (*odfValue) = _cog.cwise().clipUnder().rowwise().sum();
  
  DEBUG("ODFCOG: Finished Processing");
}

void ODFCOG::reset() {
  // Initial values
  _peaker.reset();
  _peakCoger.reset();

}
