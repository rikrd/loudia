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

#ifndef FILTER_H
#define FILTER_H

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "typedefs.h"

//using namespace std;

// import most common Eigen types
//using namespace Eigen;

class Filter {
protected:
  // Internal parameters
  Real _samplerate;
  int _channels; 
  int _length; 

  // Internal variables
  MatrixXR _a;
  MatrixXR _b;

  MatrixXR _z;
  MatrixXR _samples;

public:
  Filter(MatrixXR b, MatrixXR a, Real samplerate, int channels);

  ~Filter();

  void setup();

  void process(MatrixXR samples, MatrixXR* output);

  void reset();

  int channels() const;

  Real samplerate() const;

  int length() const;

};

#endif  /* FILTER_H */
