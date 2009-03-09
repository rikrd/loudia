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

#ifndef ODFBASE_H
#define ODFBASE_H

#include "Typedefs.h"
#include "Debug.h"

class ODFBase {
protected:
  // Internal parameters

  // Internal variables
public:
  virtual void setup() = 0;
  
  virtual void process(const MatrixXC& fft, MatrixXR* odfValue) = 0;
  
  virtual void reset() = 0;

  virtual ~ODFBase() {};
};

#endif  /* ODFBASE_H */
