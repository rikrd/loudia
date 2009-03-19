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

#ifndef IIRFILTER_H
#define IIRFILTER_H

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"
#include "FilterUtils.h"

/**
  * @class IIRFilter
  *
  * @brief Algorithm to create and apply several types of IIR filters.
  *
  * This class represents an object to create and apply several types of IIR filters.
  * Additionally the coefficients, zeros, poles and gains of the created filters
  * can be retrieved.
  *
  * 4 types of bands are implemented:
  * -# Low Pass
  * -# High Pass
  * -# Band Pass
  * -# Band Stop
  *
  * The band type can be selected using the 
  * IIRFilter::setBandType() taking as argument a BandType.
  *
  * The critical frequencies are specified using 
  * IIRFilter::setLowFrequency() and IIRFilter::setHighFrequency().
  * Note that for low pass and high pass filters which have one single critical frequency
  * only IIRFilter::setLowFrequency() has an effect.
  *
  * 4 types of filters are implemented:
  * -# Chebyshev I
  * -# Chebyshev II
  * -# Bessel
  * -# Butterworth
  *
  * The filter type can be selected using the 
  * IIRFilter::setFilterType() taking as argument a FilterType.
  *
  * The order of the filters can be specified using IIRFilter::setOrder().
  *
  * For Chebyshev I filters the pass band ripple can be specified using
  * IIRFilter::setPassRipple().  Note that this method has no effect if
  * a different type of filter is used.
  * 
  * For Chebyshev II filters the stop band attenuation is specified using
  * IIRFilter::setStopAttenuation().  Note that this method has no effect if
  * a different type of filter is used.
  *
  * @author Ricard Marxer
  *
  * @sa LowPass, HighPass, BandStop
  */
class IIRFilter {
public:
  enum FilterType {
    CHEBYSHEVI = 0,
    CHEBYSHEVII = 1,
    BUTTERWORTH = 2,
    BESSEL = 3
  };

  enum BandType {
    LOWPASS = 0,
    HIGHPASS = 1,
    BANDPASS = 2,
    BANDSTOP = 3
  };

protected:
  int _order;
  Real _lowFrequency;
  Real _highFrequency;
  Real _passRipple;
  Real _stopAttenuation;
  int _channels;
  
  Filter _filter;

  FilterType _filterType;
  BandType _bandType;
  
public:
  /**
     Constructs a band pass filter object with the given @a order, @a lowFrequency,
     @a highFrequency, @a filterType, @a ripplePass and @a attenuationStop parameters
     given.
  */
  IIRFilter(int order = 4, Real lowFrequency = 0.0, Real highFrequency = 1.0, BandType bandType = LOWPASS, FilterType filterType = CHEBYSHEVII, Real ripplePass = 0.05, Real attenuationStop = 40.0);

  void setup();
  void reset();

  void process( const MatrixXR& samples, MatrixXR* filtered );

  /**
     Return in @a a the single column matrix @a a coefficients.

     @sa b
  */
  void a( MatrixXR* a ) const;
  
  /**
     Return in @a b the single column matrix @a b coefficients.

     @sa a
  */
  void b( MatrixXR* b ) const;

  /**
     Return the order of the filter.
     The default is 4.

     @sa setOrder
  */
  int order() const;

  /**
     Specifies the @a order of the filter.
     The given @a order must be higher than 0.
     Note that orders higher than 25 are not allowed for Bessel filters.
     
     @sa order
  */
  void setOrder( int order );

  /**
     Return the low frequency of the filter.
     The default is 0.0.

     @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
  */  
  Real lowFrequency() const;  

  /**
     Specifies the low normalized @a frequency of the filter.
     The given @a frequency must be in the range of 0 to 1.
     
     @sa lowFrequency, highFrequency, setHighFrequency
  */
  void setLowFrequency( Real frequency );

  /**
     Return the stop frequency of the filter.
     The default is 1.0.

     @sa lowFrequency, setLowFrequency, setHighFrequency
  */  
  Real highFrequency() const;  

  /**
     Specifies the stop normalized @a frequency of the filter.
     The given @a frequency must be in the range of 0 to 1.
     
     @sa lowFrequency, highFrequency, setLowFrequency
  */
  void setHighFrequency( Real frequency );

  /**
     Return the filter type.
     The default is CHEBYSHEVII.
     The given @a frequency must be in the range of 0 to 1.

     @sa setFilterType, order, setOrder
  */
  FilterType filterType() const;

  /**
     Specifies the filter @a type.
     
     @sa lowFrequency, highFrequency, setLowFrequency
  */
  void setFilterType( FilterType type );

  /**
     @property IIRFilter::bandType
     @brief the type of the band of the filter
     
     By default it is LOWPASS.

     @sa filterType
  */
  BandType bandType() const;
  void setBandType( BandType type );

  /**
     @property IIRFilter::passRipple
     @brief the ripple of the pass band in dB
     
     Note that this property only has an effect if 
     the filter type used is CHEBYSHEVI.
     By default it is 0.05.

     @sa stopAttenuation
  */
  Real passRipple() const;
  void setPassRipple( Real rippleDB );

  /**
     @property IIRFilter::stopAttenuation
     @brief the attenuation of the stop band in dB
     
     Note that this property only has an effect if 
     the filter type used is CHEBYSHEVII.
     By default it is 40.0.

     @sa passRipple
  */
  Real stopAttenuation() const;
  void setStopAttenuation( Real attenuationDB );
};

#endif  /* BANDPASS_H */
