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

#include "Resample.h"

using namespace std;
using namespace Eigen;

Resample::Resample(int inSize, int outSize, Real resampleRatio, ResampleType resampleType) :
  _inSize( inSize ),
  _outSize( outSize ),
  _resampleRatio( resampleRatio ),
  _resampleType( resampleType )
{
  DEBUG("RESAMPLE: Constructor inSize: " << inSize 
        << ", outSize: " << outSize 
        << ", resampleRatio: " << resampleRatio);
  
  setup();
  
  DEBUG("RESAMPLE: Constructed");
}

Resample::~Resample(){
  DEBUG("RESAMPLE: Destroying...");

  delete [] _resampleData.data_in;
  delete [] _resampleData.data_out;
  
  DEBUG("RESAMPLE: Destroyed out");
}

void Resample::setup(){
  DEBUG("RESAMPLE: Setting up...");

  _resampleData.input_frames = _inSize;
  _resampleData.output_frames = _outSize;
  _resampleData.src_ratio = _resampleRatio;

  if ( !_resampleData.data_in ) delete [] _resampleData.data_in;
  if ( !_resampleData.data_out ) delete [] _resampleData.data_out;

  _resampleData.data_in = new float[_inSize];
  _resampleData.data_out = new float[_outSize];

  
  DEBUG("RESAMPLE: Finished set up...");
}

void Resample::process(const MatrixXR& in, MatrixXR* out){
  const int rows = in.rows();
  const int cols = in.cols();

  if ( cols != _inSize ) {
    // Throw ValueError, incorrect input size
  }

  (*out).resize(rows, _outSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer
    Eigen::Map<MatrixXR>(_resampleData.data_in, 1, _inSize) = in;
    
    // Process the data
    int error = src_simple(&_resampleData, _resampleType, 1);
    if ( error ) {
      // Throw ResampleError, src_strerror( error );
    }
    
    // Take the data from _out
    (*out).row( i ) = Eigen::Map<MatrixXR>(_resampleData.data_out, 1, _outSize);
  }
}

void Resample::reset(){
}

int Resample::inSize() const{
  return _inSize;
}

int Resample::outSize() const{
  return _outSize;
}

Real Resample::resampleRatio() const{
  return _resampleRatio;
}
