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

#ifndef BANDPASS_H
#define BANDPASS_H

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"
#include "FilterUtils.h"

/**
  * @class BandPass
  *
  * @brief Algorithm to create and apply band pass filters.
  *
  * This class represents an object to create and apply band pass filters.
  * Additionally the coefficients, zeros, poles and gains of the created filters
  * can be retrieved.
  *
  * 4 types of filters are implemented:
  * -# Chebyshev I
  * -# Chebyshev II
  * -# Bessel
  * -# Butterworth
  *
  * For Chebyshev I filters 
  *
  * Note that the number of rows of the starts matrix and the size of the vector of weights must
  * be the same, and this will be the number of bands.
  *
  * @author Ricard Marxer
  *
  * @sa MelBands
  */
class BandPass {
protected:
  int _order;
  Real _startFrequency;
  Real _stopFrequency;
  Real _passRipple;
  Real _stopAttenuation;
  int _channels;
  
  Filter _filter;

  FilterType _filterType;

public:
  BandPass(int order = 4, Real startFrequency = 0.2, Real stopFrequency = 0.4, FilterType filterType = CHEBYSHEVII, Real ripplePass = 0.05, Real attenuationStop = 40.0);

  void setup();
  void reset();

  void process( const MatrixXR& samples, MatrixXR* filtered );

  void a( MatrixXR* a ) const;
  void b( MatrixXR* b ) const;
  
  int order() const;
  void setOrder( int order );

  Real startFrequency() const;  
  void setStartFrequency( Real frequency );

  Real stopFrequency() const;  
  void setStopFrequency( Real frequency );

  FilterType filterType() const;
  void setFilterType( FilterType type );

  Real passRipple() const;
  void setPassRipple( Real rippleDB );

  Real stopAttenuation() const;
  void setStopAttenuation( Real attenuationDB );
};

#endif  /* BANDPASS_H */
