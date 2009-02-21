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

#ifndef BANDS_H
#define BANDS_H

#include <Eigen/StdVector>
#include <vector>

#include "Typedefs.h"
#include "Debug.h"

class Bands {
protected:
  // Internal parameters
  MatrixXI _starts;
  std::vector<MatrixXR> _weights;

  // Internal variables

public:
  Bands();

  Bands(MatrixXI starts, std::vector<MatrixXR> weights);

  ~Bands();

  void setup();

  void process(const MatrixXR&  spectrum, MatrixXR* bands);

  void reset();

  std::vector<MatrixXR> weights() const;

  void bandWeights(int band, MatrixXR* bandWeights) const;

  void starts(MatrixXI* result) const;

  int bands() const;

  void setStartsWeights(const MatrixXI& starts, std::vector<MatrixXR> weights);
};

#endif  /* BANDS_H */