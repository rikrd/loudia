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

#include "SpectralODF.h"
#include "SpectralODFComplex.h"
#include "SpectralODFPhase.h"
#include "SpectralODFMKL.h"
#include "SpectralODFFlux.h"
#include "SpectralODFHFC.h"
#include "SpectralODFCOG.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralODF::SpectralODF(int fftSize, ODFMethod odfMethod) :
  SpectralODFBase(),
  _odf( 0 )
{
  setFftSize( fftSize, false );
  setOdfMethod( odfMethod, false );

  setup();
}

SpectralODF::~SpectralODF() {
  delete _odf;
  _odf = 0;
}

void SpectralODF::setup() {

  delete _odf;
  _odf = 0;

  switch( _odfMethod ) {

  case FLUX:
    _odf = new SpectralODFFlux(_fftSize);
    break;

  case PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize);
    break;

  case WEIGHTED_PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize, true);
    break;

  case NORM_WEIGHTED_PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize, true, true);
    break;

  case MODIFIED_KULLBACK_LIEBLER:
    _odf = new SpectralODFMKL(_fftSize);
    break;

  case COMPLEX_DOMAIN:
    _odf = new SpectralODFComplex(_fftSize);
    break;

  case RECTIFIED_COMPLEX_DOMAIN:
    _odf = new SpectralODFComplex(_fftSize, true);
    break;

  case HIGH_FREQUENCY_CONTENT:
    _odf = new SpectralODFHFC(_fftSize);
    break;

  case CENTER_OF_GRAVITY:
    _odf = new SpectralODFCOG(_fftSize);
    break;

  }

  _odf->setup();
}

void SpectralODF::process(const MatrixXC& fft, MatrixXR* odfValue) {
  if (!_odf) return;
  
  _odf->process(fft, odfValue);
}

void SpectralODF::reset() {
  if (!_odf) return;
  
  _odf->reset();
}

SpectralODF::ODFMethod SpectralODF::odfMethod() const{
  return _odfMethod;
}

void SpectralODF::setOdfMethod( ODFMethod method, bool callSetup ) {
  _odfMethod = method;
  if ( callSetup ) setup();
}
