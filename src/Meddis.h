/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#ifndef MEDDIS_H
#define MEDDIS_H

#include "Typedefs.h"
#include "Debug.h"

class Meddis {
protected:
  // Internal fixed constants
  static const Real M;
  static const Real A;
  static const Real B;
  static const Real g;
  static const Real y;
  static const Real l;
  static const Real r;
  static const Real x;
  static const Real h;

  // Internal parameters
  Real _samplerate;
  int _channels;
  bool _substractSpont;
  
  Real dt;
  Real gdt;
  Real ydt;
  Real ldt;
  Real rdt;
  Real xdt;

  // Internal variables
  MatrixXR kt;
  MatrixXR spont;

  MatrixXR c;
  MatrixXR q;
  MatrixXR w;

public:
  Meddis(Real samplerate, int channels, bool substractSpont = true);

  ~Meddis();

  void setup();

  void process(const MatrixXR& data, MatrixXR* output);

  int channels() const;

  Real samplerate() const;

  void reset();
};

#endif  /* MEDDIS_H */
