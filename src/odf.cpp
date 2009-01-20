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

#include "odf.h"
#include "odfcomplex.h"

#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODF::ODF(int fftLength, ODFType odfType) :
  _odfType(odfType),
  _fftLength(fftLength)
{
  switch(_odfType) {

  case SPECTRAL_FLUX:
  case PHASE_DEVIATION:
  case WEIGHTED_PHASE_DEVIATION:
  case NORM_WEIGHTED_PHASE_DEVIATION:
  case MODIFIED_KULLBACK_LIEBLER:
    // Throw ImplementationError, ODF type not implemented yet
    break;

  case COMPLEX_DOMAIN:
    _odf = new ODFComplex(_fftLength);
    break;

  case RECTIFIED_COMPLEX_DOMAIN:
    break;
  }
  
}

ODF::~ODF() {
  delete _odf;  
}

void ODF::setup() {
  _odf->setup();
}

void ODF::process(const MatrixXC& fft, MatrixXR* odfValue) {
  _odf->process(fft, odfValue);
}

void ODF::reset() {
  _odf->reset();
}
