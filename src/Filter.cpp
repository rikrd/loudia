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

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"

using namespace std;
using namespace Eigen;

Filter::Filter(int channels) : 
  _channels( channels )
{

  setB(MatrixXR::Ones(1, _channels));

  setA(MatrixXR::Ones(1, _channels));

  setup();
}


Filter::Filter(const MatrixXR& b,
               const MatrixXR& a,
               int channels) :
  _channels( channels ),
  _ina( a ),
  _inb( b )
{

  DEBUG("FILTER: Constructor channels:" << channels);
  DEBUG("FILTER: Constructor b:" << b.transpose() << ", a:" << a.transpose());
    
  setup();
}

Filter::~Filter() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void Filter::setup(){
  // Prepare the buffers
  DEBUG("FILTER: Setting up...");

  _length = max(_inb.rows(), _ina.rows());
  
  _z.resize(_length, _channels);

  _samples.resize(1, _channels);

  // Normalize by the first value value the denominator coefficients
  // since a[0] must be 1.0
  // TODO: throw an exception when a[0] == 0
  DEBUG("FILTER: Initializing 'a' coeffs");
  _a = MatrixXR::Zero(_length, _channels);

  // Check that it has one column or as many as channels
  if ((_ina.cols() != 1) && (_ina.cols() != _channels)) {
    // TODO: Throw an exception
    DEBUG("FILTER: Error in shape of 'a' coeffs. _ina.cols():" << _ina.cols() << ", _channels:" << _channels);
    return;
  }

  // Set the coefficients
  if (_ina.cols() == 1) {
    // If only one column has been defined repeat it for all columns
    for (int i=0; i < _a.cols(); i++) {
      _a.block(0, i, _ina.rows(), 1) = _ina.col(0);
    }
  }else{
    // Else set it directly
    _a.block(0, 0, _ina.rows(), _channels) = _ina;
  }

  for(int i = 0; i < _a.rows(); i++){
    _a.row(i) = _a.row(i).cwise() / _ina.row(0);
  }

  DEBUG("FILTER: Setting the 'a' coefficients.");
  DEBUG("FILTER: 'a' coefficients: " << _a.transpose());
  
  _b = MatrixXR::Zero(_length, _channels);

  // Check that it has one column or as many as channels
  if ((_inb.cols() != 1) && (_inb.cols() != _channels)) {
    // TODO: Throw an exception
    DEBUG("FILTER: Error in shape of 'b' coeffs. b.cols():" << _inb.cols() << ", _channels:" << _channels);
    return;
  }

  // Set the coefficients
  if (_inb.cols() == 1) {
    // If only one column has been defined repeat it for all columns
    for (int i=0; i < _b.cols(); i++) {
      _b.block(0, i, _inb.rows(), 1) = _inb.col(0);
    }
  }else{
    // Else set it directly
    _b.block(0, 0, _inb.rows(), _channels) = _inb;
  }

  for(int i = 0; i < _b.rows(); i++){
    _b.row(i) = _b.row(i).cwise() / _ina.row(0);
  }  
  
  DEBUG("FILTER: Setting the 'b' coefficients.");
  DEBUG("FILTER: 'b' coefficients: " << _b.transpose());  
  
  reset();

  DEBUG("FILTER: Finished set up...");
}


void Filter::process(const MatrixXR& samples, MatrixXR* output){
  // Process will be called with a matrix where columns will be channels
  // and rows will be the time axis

  //DEBUG("FILTER: Entering process...");
  _samples.resize(samples.rows(), _channels);
  (*output).resize(samples.rows(), _channels);
  //DEBUG("FILTER: After resize...");  
  

  // Check that it has one column or as many as channels
  if ((samples.cols() != 1) && (samples.cols() != _channels)) {
    // TODO: Throw an exception
    DEBUG("FILTER: Error in shape of 'samples'. samples.cols():" << samples.cols() << ", _channels:" << _channels);
    return;
  }

  // Set the input
  if (samples.cols() == 1) {
    // If only one column has been defined repeat it for all columns
    for (int i=0; i < _samples.cols(); i++) {
      _samples.block(0, i, samples.rows(), 1) = samples.col(0);
    }
  }else{
    // Else set it directly
    _samples.block(0, 0, samples.rows(), _channels) = samples;
  }
  //DEBUG("FILTER: _a: " << _a);

  //DEBUG("FILTER: After setting coeffs...");    
  for ( int i = 0; i < _samples.rows(); i++ ) {
    if ( _length > 1 ) {
      //DEBUG("output.rows(): " << (*output).rows());
      //DEBUG("output.cols(): " << (*output).cols());
      //DEBUG("_z.rows(): " << _z.rows());
      //DEBUG("_z.cols(): " << _z.cols());
      //DEBUG("_b.rows(): " << _b.rows());
      //DEBUG("_b.cols(): " << _b.cols());
      //DEBUG("_a.rows(): " << _a.rows());
      //DEBUG("_a.cols(): " << _a.cols());
      (*output).row( i ) = _z.row( 0 ) + (_b.row( 0 ).cwise() * _samples.row( i ));
      
      //DEBUG("FILTER: After setting output..., output: " << (*output).row( i ));
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
  //DEBUG("FILTER: output: " << (*output));
  //DEBUG("FILTER: After processing...");
}


void Filter::setA(const MatrixXR& a){
  _ina = a;
}

void Filter::setB(const MatrixXR& b){
  _inb = b;
}

void Filter::a(MatrixXR* a){
  (*a) = _a;
}

void Filter::b(MatrixXR* b){
  (*b) = _b;
}


void Filter::reset(){
  // Initial values
  _z = MatrixXR::Zero(_length, _channels);
}

int Filter::channels() const {
  return _channels;
}

int Filter::length() const {
  return _length;
}
