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

#include "VoiceActivityDetection.h"

using namespace std;
using namespace Eigen;

VoiceActivityDetection::VoiceActivityDetection(int lowBand, int highBand,
                                               Real sampleRate,
                                               int fftSize,
                                               int memorySize) : 
  _lowBand(lowBand),
  _highBand(highBand),
  _sampleRate(sampleRate),
  _fftSize(fftSize),
  _memorySize(memorySize),
  _barkBands(lowBand, highBand, sampleRate, fftSize)
{

  LOUDIA_DEBUG("VoiceActivityDetection: Constructor");

  setup();
  
  LOUDIA_DEBUG("VoiceActivityDetection: Constructed");
}

VoiceActivityDetection::~VoiceActivityDetection(){
  LOUDIA_DEBUG("VoiceActivityDetection: Destroying...");
  LOUDIA_DEBUG("VoiceActivityDetection: Destroyed out");
}

void VoiceActivityDetection::setup(){
  LOUDIA_DEBUG("VoiceActivityDetection: Setting up...");

  _memory = MatrixXR::Zero(_memorySize, _highBand - _lowBand + 1);
  _currentMemoryPos = 0;
  
  _barkBands.setSampleRate(_sampleRate, false);
  _barkBands.setLowBand(_lowBand, false);
  _barkBands.setHighBand(_highBand, false);
  _barkBands.setup();
  
  LOUDIA_DEBUG("VoiceActivityDetection: Finished set up...");
}

void VoiceActivityDetection::process(const MatrixXR& frames, MatrixXR* vad){
  const int rows = frames.rows();
  
  vad->resize(rows, 1);

  for (int i=0; i < rows; i++){
    // compute barkbands
    _barkBands.process(frames.row(0), &_bands);

    // copy frame into memory
    _memory.row(_currentMemoryPos) = _bands.row(0);

    _currentMemoryPos = (_currentMemoryPos + 1) % _memorySize;

    // compute the VAD
    RowXR LTSE = _memory.colwise().maxCoeff();
    RowXR noise = _memory.colwise().sum() / _memorySize;

    (*vad)(i,0) = log10((LTSE.array().square() / noise.array().square()).sum());
  }
}

void VoiceActivityDetection::reset(){
  _barkBands.reset();
}
