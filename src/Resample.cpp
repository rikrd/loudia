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

Resample::Resample(int inputSize, int outputSize, Real resamplingRatio, ResamplingMethod resamplingMethod)
{
  DEBUG("RESAMPLE: Constructor inputSize: " << inputSize 
        << ", outputSize: " << outputSize 
        << ", resamplingRatio: " << resamplingRatio);

  _resampleData.data_in = NULL;
  _resampleData.data_out = NULL;
  
  setInputSize( inputSize, false );
  setOutputSize( outputSize, false );
  setResamplingRatio( resamplingRatio, false );
  setResamplingMethod( resamplingMethod, false );

  setup();
  
  DEBUG("RESAMPLE: Constructed");
}

Resample::~Resample(){
  DEBUG("RESAMPLE: Destroying...");

  if ( _resampleData.data_in ) delete [] _resampleData.data_in;
  if ( _resampleData.data_out ) delete [] _resampleData.data_out;
  
  DEBUG("RESAMPLE: Destroyed out");
}

void Resample::setup(){
  DEBUG("RESAMPLE: Setting up...");

  _resampleData.input_frames = _inputSize;
  _resampleData.output_frames = _outputSize;
  _resampleData.src_ratio = _resamplingRatio;

  if ( _resampleData.data_in ) delete [] _resampleData.data_in;
  if ( _resampleData.data_out ) delete [] _resampleData.data_out;

  _resampleData.data_in = new float[_inputSize];
  _resampleData.data_out = new float[_outputSize];

  
  DEBUG("RESAMPLE: Finished set up...");
}

void Resample::process(const MatrixXR& in, MatrixXR* out){
  const int rows = in.rows();
  const int cols = in.cols();

  if ( cols != _inputSize ) {
    // Throw ValueError, incorrect input size
  }

  (*out).resize(rows, _outputSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer
    Eigen::Map<MatrixXR>(_resampleData.data_in, 1, _inputSize) = in;
    
    // Process the data
    int error = src_simple(&_resampleData, _resamplingMethod, 1);
    if ( error ) {
      // Throw ResampleError, src_strerror( error );
    }
    
    // Take the data from _out
    (*out).row( i ) = Eigen::Map<MatrixXR>(_resampleData.data_out, 1, _outputSize);
  }
}

void Resample::reset(){
}

int Resample::inputSize() const{
  return _inputSize;
}

void Resample::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}

int Resample::outputSize() const{
  return _outputSize;
}

void Resample::setOutputSize( int size, bool callSetup ) {
  _outputSize = size;
  if ( callSetup ) setup();
}

Real Resample::resamplingRatio() const{
  return _resamplingRatio;
}

void Resample::setResamplingRatio( Real ratio, bool callSetup ) {
  _resamplingRatio = ratio;
  if ( callSetup ) setup();
}

Resample::ResamplingMethod Resample::resamplingMethod() const{
  return _resamplingMethod;
}

void Resample::setResamplingMethod( ResamplingMethod method, bool callSetup ) {
  _resamplingMethod = method;
  if ( callSetup ) setup();
}
