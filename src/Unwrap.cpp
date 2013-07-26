/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

#include "Typedefs.h"
#include "Debug.h"

#include "Unwrap.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

Unwrap::Unwrap()
{
  LOUDIA_DEBUG("UNWRAP: Construction");

  setup();
}

Unwrap::~Unwrap(){}

void Unwrap::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("UNWRAP: Setting up...");
  
  reset();

  LOUDIA_DEBUG("UNWRAP: Finished setup.");
}

void Unwrap::process(const MatrixXR& input, MatrixXR* unwrapped){
  const int rows = input.rows();
  const int cols = input.cols();

  (*unwrapped).resize(rows, cols);
  
  if(input.rows() <= 1){
    (*unwrapped) = input;
  }

  _diff.resize(rows, cols);
  _upsteps.resize(rows, cols);
  _downsteps.resize(rows, cols);
  _shift.resize(rows, cols);  

  _diff << MatrixXR::Zero(1, cols), input.block(0, 0, rows-1, cols) - input.block(1, 0, rows-1, cols);
  
  _upsteps = (_diff.array() > M_PI).cast<Real>();
  _downsteps = (_diff.array() < -M_PI).cast<Real>();

  rowCumsum(&_upsteps);
  rowCumsum(&_downsteps);

  _shift =  _upsteps - _downsteps;

  (*unwrapped) = input + (2.0 * M_PI * _shift);
}

void Unwrap::reset(){
  // Initial values
}
