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

#include "FrameCutter.h"

#include <iostream>
using namespace std;

FrameCutter::FrameCutter( const int maxInputSize, const int frameSize, const int hopSize, const int firstSamplePosition, const Real defaultValue ) :
  _defaultValue(defaultValue)
{
  LOUDIA_DEBUG("FRAMECUTTER: Constructing...");

  setMaxInputSize(maxInputSize, false);

  LOUDIA_DEBUG("FRAMECUTTER: Set the maximum input size...");

  setFrameSize(frameSize, false);

  LOUDIA_DEBUG("FRAMECUTTER: Set the frame size...");

  setHopSize(hopSize, false);

  LOUDIA_DEBUG("FRAMECUTTER: Set the hop size...");

  setFirstSamplePosition(firstSamplePosition, false);

  LOUDIA_DEBUG("FRAMECUTTER: Set first sample position...");

  setup();
}

FrameCutter::~FrameCutter(){
  
}

void FrameCutter::setup(){  
  const int bufferSize = 2 * _frameSize;

  // Create the stream buffer
  // must be at least twice the size than the frame size
  // in order to keep data aligned
  _buffer.resize(bufferSize);
  _buffer.setConstant(_defaultValue);
  
  if ((_firstSamplePosition < 0) || (_firstSamplePosition > _frameSize - 1)) {
    LOUDIA_ERROR("FRAMECUTTER: The first sample position must be set between 0 and frameSize-1.");
  }
  
  _indexWriter = _firstSamplePosition;
  _availableToWrite = _frameSize - _firstSamplePosition;
  _availableToRead = _firstSamplePosition;

  // TODO: check if this is the right way to know the maxFrameCount
  _maxFrameCount = maxFrameCount();

  _row = VectorXR::Zero(_frameSize);
}

void FrameCutter::process(const MatrixXR& stream, MatrixXR *frames, int *produced){
  if (stream.cols() != 1) {
    LOUDIA_ERROR("FRAMECUTTER: This algorithm only accepts single channel streams.");
  }

  frames->resize(_maxFrameCount, _frameSize);
  
  int currenthopSize = hopSize();
  
  int leftSize = stream.rows();
  int inputIndex = 0;
  int inputSize;
  int framesIndex = 0;
  
  while (leftSize > 0) {
    int consumed = 1;    
    while (consumed > 0) {
      inputSize = min(_frameSize, leftSize);
      if (inputSize <= 0) break;
    
      consumed = write(stream.col(0).segment(inputIndex, inputSize));
      leftSize -= consumed;
      inputIndex += consumed;
    }  
    
    while (read(&_row, currenthopSize) > 0) {
      // TODO: try to avoid the copy to frames->row(framesIndex)
      // Maybe by passing directly frames->row(framesIndex) as a const ref
      // and const casting in read()
      frames->row(framesIndex) = _row;
      framesIndex += 1;
    }
  }
  
  (*produced) = framesIndex;
  return;
}

int FrameCutter::read(VectorXR* frame, int release){
  if ( frame->size() > _availableToRead) return 0;
  
  const int indexReader = ((_indexWriter - _availableToRead) + _frameSize) % _frameSize;
  (*frame) = _buffer.segment(indexReader, frame->size());
  
  // Advance reader index  
  _availableToRead -= release;
  _availableToWrite += release;
  
  return frame->size();
}

int FrameCutter::write(const VectorXR& stream){
  int consumed = min(min(_availableToWrite, _frameSize - _indexWriter), (int)stream.size());
  if ( consumed <= 0 ) return 0;

  _buffer.segment(_indexWriter, consumed) = stream.segment(0, consumed);
  _buffer.segment(_indexWriter + _frameSize, consumed) = stream.segment(0, consumed);
  
  // Advance writer index
  _indexWriter = (_indexWriter + consumed) % _frameSize;
  _availableToWrite -= consumed;
  _availableToRead += consumed;
  
  return consumed;
}

void FrameCutter::reset(){
  _indexWriter = 0;
  _availableToWrite = _frameSize;
  _availableToRead = 0;
}

int FrameCutter::maxFrameCount() const {
  return _maxInputSize / hopSize() + 1; 
}
