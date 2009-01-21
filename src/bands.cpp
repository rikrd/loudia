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

#include "typedefs.h"
#include "debug.h"

#include "bands.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

Bands::Bands() : _starts(1, 1){ 
  _weights.push_back(MatrixXR::Constant(1, 1, 0.0));

  setup();
}


Bands::Bands(MatrixXi starts, vector<MatrixXR> weights) {
  DEBUG("BANDS: Constructor starts: " << starts);

  if ( starts.rows() != (int)weights.size() ) {
    // Throw an exception
  }

  if ( starts.cols() != 1 ) {
    // Throw an exception
  }


  for (int i = 0; i < (int)weights.size(); i++){
    if ( weights[i].cols() != 1 ) {
      // Throw an exception
    }
  }

  _starts = starts;
  _weights = weights;

  setup();
}

Bands::~Bands() {}

void Bands::setup(){
  // Prepare the buffers
  DEBUG("BANDS: Setting up...");

  reset();
  DEBUG("BANDS: Finished set up...");
}


void Bands::process(const MatrixXR& spectrum, MatrixXR* bands){

  (*bands).resize(spectrum.rows(), _starts.rows());

  for (int j = 0; j < spectrum.rows(); j++){
    for (int i = 0; i < _starts.rows(); i++ ) {
      (*bands)(j, i) = spectrum.block(j, _starts(i, 0), 1, _weights[i].rows()).row(0).dot(_weights[i].col(0));
    }
  }
}

void Bands::reset(){
  // Initial values
}

vector<MatrixXR> Bands::weights() const {
  return _weights;
}

void Bands::bandWeights(int band, MatrixXR* bandWeights) const {
  (*bandWeights) =  _weights[ band ];
}

void Bands::starts(MatrixXI* result) const {
  (*result) = _starts;
}

void Bands::setStartsWeights(const MatrixXI& starts, std::vector<MatrixXR> weights) {
  _weights = weights;
  _starts = starts;
}

int Bands::bands() const {
  return _starts.rows();
}
