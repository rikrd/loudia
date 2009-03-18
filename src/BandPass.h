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
  * The filter type can be selected using the BandPass::setFilterType() taking as
  * argument a FilterType.
  *
  * The order of the filters can be specified using BandPass::setOrder().
  * The critical frequencies are specified using 
  * BandPass::setStartFrequency() and BandPass::setStopFrequency(). 
  *
  * For Chebyshev I filters the pass band ripple can be specified using
  * BandPass::setPassRipple().  Note that this method has no effect if
  * a different type of filter is used.
  * 
  * For Chebyshev II filters the stop band attenuation is specified using
  * BandPass::setStopAttenuation().  Note that this method has no effect if
  * a different type of filter is used.
  *
  * @author Ricard Marxer
  *
  * @sa LowPass, HighPass, BandStop
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
  /**
     Constructs a band pass filter object with the given @a order, @a startFrequency,
     @a stopFrequency, @a filterType, @a ripplePass and @a attenuationStop parameters
     given.
  */
  BandPass(int order = 4, Real startFrequency = 0.2, Real stopFrequency = 0.4, FilterType filterType = CHEBYSHEVII, Real ripplePass = 0.05, Real attenuationStop = 40.0);

  void setup();
  void reset();

  void process( const MatrixXR& samples, MatrixXR* filtered );

  /**
     Return in @a a the single column matrix @a a coefficients.
  */
  void a( MatrixXR* a ) const;
  
  /**
     Return in @a b the single column matrix @a b coefficients.
  */
  void b( MatrixXR* b ) const;

  /**
     Return the order of the filter.
     The default is 4.
  */  
  int order() const;
  void setOrder( int order );

  /**
     Return the start frequency of the filter.
     The default is 0.2.
  */  
  Real startFrequency() const;  
  void setStartFrequency( Real frequency );

  /**
     Return the stop frequency of the filter.
     The default is 0.4.
  */  
  Real stopFrequency() const;  
  void setStopFrequency( Real frequency );

  /**
     Return the filter type.
     The default is CHEBYSHEVII.
  */
  FilterType filterType() const;
  void setFilterType( FilterType type );

  /**
     Return the pass ripple in dB.
     The default is 0.05.
  */
  Real passRipple() const;
  void setPassRipple( Real rippleDB );

  /**
     Return the stop attenuation in dB.
     The default is 40.0.
  */
  Real stopAttenuation() const;
  void setStopAttenuation( Real attenuationDB );
};

#endif  /* BANDPASS_H */
