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
#include "filter.h"

#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

Filter::Filter(MatrixXR b,
               MatrixXR a,
               Real samplerate,
               int channels) {

  DEBUG("FILTER: Constructor {samplerate:" << samplerate << ", channels:" << channels << "}")
  
  _channels = channels;
  _length = max(b.rows(), a.rows());
  _samplerate = samplerate;

  _samples(1, _channels);
  //cout << "Building..."<<endl;

  // Normalize by the first value value the denominator coefficients
  // since a[0] must be 1.0
  // TODO: throw an exception when a[0] == 0
  _a = MatrixXR::Zero(_length, _channels);
  _a.block(0, 0, a.rows(), _channels) = a;

  for(int i = 0; i < _a.rows(); i++){
    _a.row(i) = _a.row(i).cwise() / _a.row(0);
  }

  //cout << "a coeffs:" << _a << endl;

  //cout << "built _a..."<<endl;  
  _b = MatrixXR::Zero(_length, _channels);
  _b.block(0, 0, b.rows(), _channels) = b;

  for(int i = 0; i < _b.rows(); i++){
    _b.row(i) = _b.row(i).cwise() / _a.row(0);
  }  
  
  //cout << "b coeffs:" << _b << endl;
  //cout << "built _b..."<<endl;  
}

Filter::~Filter() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void Filter::setup(){
  // Prepare the buffers
  DEBUG("FILTER: Setting up...")

  reset();
}


void Filter::process(MatrixXR samples, MatrixXR* output){
  // Process will be called with a matrix where columns will be channels
  // and rows will be the time axis
  
  _samples.resize(samples.rows(), _channels);
  
  // If there is one single input channel then repeat it 
  if (samples.cols() != _channels){
    for (int i= 0; i < _samples.cols(); i++) {
      _samples.col(i) = samples.col(0);
      
    }
  } else {
    _samples = samples;
  }

  for ( int i = 0; i < samples.rows(); i++ ) {

    if ( _length > 1 ) {
      (*output).row( i ) = _z.row( 0 ) + (_b.row( 0 ).cwise() * _samples.row( i ));
      
      // Fill in middle delays
      for ( int j = 0; j < _length - 1; j++ ) {
        _z.row( j ) = _z.row( j + 1 ) + (_samples.row( i ).cwise() * _b.row( j + 1 )) - ((*output).row( i ).cwise() * _a.row( j + 1 ));
      }

      // Calculate the last delay
      _z.row( _length - 2 ) = (_samples.row( i ).cwise() * _b.row( _length - 1 )) - ((*output).row( i ).cwise() * _a.row( _length - 1 ));
    } else {
      (*output).row( i ) = _samples.row( i ) * _b.row( 0 );
    }
  }  
}

void Filter::reset(){
  // Initial values
  _z = MatrixXR::Zero(_length, _channels);
}

int Filter::channels() const {
  return _channels;
}

Real Filter::samplerate() const {
  return _samplerate;
}

int Filter::length() const {
  return _length;
}
