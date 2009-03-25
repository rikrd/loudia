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

#ifndef SPECTRALODF_H
#define SPECTRALODF_H

#include "Typedefs.h"
#include "Debug.h"

#include "SpectralODFBase.h"

class SpectralODF : SpectralODFBase {
public:
  enum SpectralODFType {
    FLUX = 0,
    HIGH_FREQUENCY_CONTENT = 1,
    PHASE_DEVIATION = 2,
    WEIGHTED_PHASE_DEVIATION = 3,
    NORM_WEIGHTED_PHASE_DEVIATION = 4,
    MODIFIED_KULLBACK_LIEBLER = 5,
    COMPLEX_DOMAIN = 6,
    RECTIFIED_COMPLEX_DOMAIN = 7,
    CENTER_OF_GRAVITY = 8
  };

protected:
  // Internal parameters
  SpectralODFType _odfType;
  
  // Internal variables
  SpectralODFBase* _odf;

public:
  SpectralODF(int fftSize, SpectralODFType odfType = COMPLEX_DOMAIN);
  
  ~SpectralODF();

  void setup();
  void reset();

  void process(const MatrixXC& fft, MatrixXR* odfValue);


};

#endif  /* SPECTRALODF_H */
