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

#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include "filter.h"

#include "typedefs.h"

class Chebyshev {
protected:
  Real _rippleDB;
  Real _samplerate;
  int _order;
  int _channels;

  Filter _filter;

public:
  Chebyshev(int channels, int order, Real rippleDB, Real samplerate);

  void setup();

  void process(MatrixXR samples, MatrixXR* filtered);

  void roots(MatrixXR poly, MatrixXC* roots);
  
  void reset();
};

#endif  /* CHEBYSHEV_H */
