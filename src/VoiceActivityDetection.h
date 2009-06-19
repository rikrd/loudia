/*                                                         
** Copyright (C) 2008, 2009 Nicolas Wack <wackou@gmail.com>
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

#ifndef VOICEACTIVITYDETECTION_H
#define VOICEACTIVITYDETECTION_H

#include "Typedefs.h"
#include "Debug.h"


/**
 * WARNING: the process() method is NOT re-entrant.
 */
class VoiceActivityDetection {

protected:
  Real _sampleRate;
  int _lowBand, _highBand;
  MatrixXR _memory;
  int _currentMemoryPos;
  int _memorySize;

  int _halfSize;
    

public:
  VoiceActivityDetection(Real sampleRate = 44100,
                         int lowBand = 4, int highBand = 16,
                         int memorySize = 12);
  ~VoiceActivityDetection();
  
  void process(const MatrixXR& frames, MatrixXR* vad);

  void setup();
  void reset();


  // Parameters
  Real sampleRate() const { return _sampleRate; }
  void setSampleRate(Real sampleRate) { _sampleRate = sampleRate; }

  Real lowBand() const { return _lowBand; }
  void setLowBand(Real lowBand) { _lowBand = lowBand; }

  Real highBand() const { return _highBand; }
  void setHighBand(Real highBand) { _highBand = highBand; }

  Real memorySize() const { return _memorySize; }
  void setMemorySize(Real memorySize) { _memorySize = memorySize; }

};

#endif  // VOICEACTIVITYDETECTION_H
