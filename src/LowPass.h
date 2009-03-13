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

#ifndef LOWPASS_H
#define LOWPASS_H

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"
#include "FiltersUtils.h"

class LowPass {
protected:
  int _order;
  Real _freq;
  Real _rippleDB;
  int _channels;
  
  Filter _filter;

  FilterType _filterType;

public:
  LowPass(int order, Real freq, Real rippleDB, int channels = 1, FilterType filterType = CHEBYSHEVI);

  void setup();

  void process(MatrixXR samples, MatrixXR* filtered);

  void a(MatrixXR* a);
  void b(MatrixXR* b);
  
  void reset();
};

#endif  /* CHEBYSHEV_H */
