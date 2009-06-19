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

VoiceActivityDetection::VoiceActivityDetection(Real sampleRate,
                                               int lowBand, int highBand,
                                               int memorySize)
  : _sampleRate(sampleRate), _lowBand(lowBand), _highBand(highBand), _memorySize(memorySize) {

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
  
  LOUDIA_DEBUG("VoiceActivityDetection: Finished set up...");
}

void VoiceActivityDetection::process(const MatrixXR& frames, MatrixXR* vad){
  const int cols = frames.cols();
  const int rows = frames.rows();
  
  //cout << "vad resize " << rows << endl;

  vad->resize(rows, 1);

  RowXR startFreqs(26);
  startFreqs << 0, 20, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500;
  
  int nbands = _highBand - _lowBand + 1;
  RowXI startBins = ((startFreqs.segment(_lowBand, nbands+1) / _sampleRate) * (cols-1)*2).cast<int>();

  RowXR bbands(nbands);

  //cout << "framesize " << rows << " - " << cols << endl;

  for (int i=0; i < rows; i++){
    //cout << "loop" << i << endl;
    // compute barkbands
    for (int b=0; b<nbands; b++) {
      //cout << " - " << startBins[b] << " " << startBins[b+1] << flush;
      bbands[b] = frames.block(i, startBins[b], 1, startBins[b+1]-startBins[b]).cwise().square().sum();
    }
    //cout << bbands << endl;
    
    // copy frame into memory
    _memory.row(_currentMemoryPos) = bbands;

    _currentMemoryPos = (_currentMemoryPos + 1) % _memorySize;

    // compute the VAD
    RowXR LTSE = _memory.colwise().maxCoeff();
    RowXR noise = _memory.colwise().sum() / _memorySize;

    (*vad)(i,0) = log10((LTSE.cwise().square().cwise() / noise.cwise().square()).sum());
  }
}

void VoiceActivityDetection::reset(){
}
