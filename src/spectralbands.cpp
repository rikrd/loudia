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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>

#include "spectralbands.h"
#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

SpectralBands::SpectralBands() : _starts(1, 1){ }


SpectralBands::SpectralBands(MatrixXi starts, vector<MatrixXR> weights) {
  if ( starts.rows() != weights.size() ) {
    // Throw an exception
  }

  if ( starts.cols() != 1 ) {
    // Throw an exception
  }


  for (int i = 0; i < weights.size(); i++){
    if ( weights[i].cols() != 1 ) {
      // Throw an exception
    }
  }

  _starts = starts;
  _weights = weights;
}

SpectralBands::~SpectralBands() {}

void SpectralBands::setup(){
  // Prepare the buffers
  reset();
}


void SpectralBands::process(MatrixXR spectrum, MatrixXR* bands){
  (*bands).resize(spectrum.rows(), _starts.rows());

  for (int j = 0; j < spectrum.rows(); j++){
    //DEBUG("SPECTRALBANDS: Process spectrum.cols(): " << spectrum.cols() << ", spectrum.rows():" << spectrum.rows() << "");

    for (int i = 0; i < _starts.rows(); i++ ) {
      //DEBUG("SPECTRALBANDS: Process _weight[" << i <<"]: [" << _weights[i].transpose() << "]");
      //DEBUG("SPECTRALBANDS: Process _weight[" << i <<"].rows(): " << _weights[i].rows() << "");
      //DEBUG("SPECTRALBANDS: Process _starts(" << i <<", 0): " << _starts(i, 0) << "");
      (*bands)(j, i) = spectrum.block(j, _starts(i, 0), 1, _weights[i].rows()).row(0).dot(_weights[i].col(0));
    }
  }
}

void SpectralBands::reset(){
  // Initial values
}

vector<MatrixXR> SpectralBands::weights() const {
  return _weights;
}

MatrixXi SpectralBands::starts() const {
  return _starts;
}
