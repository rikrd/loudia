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

#include "BarkBands.h"

/**
 * WARNING: the process() method is NOT re-entrant.
 */
class VoiceActivityDetection {

protected:
  int _lowBand, _highBand;
  Real _sampleRate;
  int _fftSize;
  MatrixXR _memory;
  int _currentMemoryPos;
  int _memorySize;

  BarkBands _barkBands;
  MatrixXR _bands;

  int _halfSize;
    

public:
  VoiceActivityDetection(int lowBand = 4, int highBand = 16,
                         Real sampleRate = 44100,
                         int fftSize = 1024,
                         int memorySize = 12);
  ~VoiceActivityDetection();
  
  void process(const MatrixXR& frames, MatrixXR* vad);

  void setup();
  void reset();


  // Parameters
  int lowBand() const { return _lowBand; }
  void setLowBand(int lowBand, bool callSetup = true ) { _lowBand = lowBand; if ( callSetup ) setup(); }

  int highBand() const { return _highBand; }
  void setHighBand(int highBand, bool callSetup = true ) { _highBand = highBand; if ( callSetup ) setup(); }

  Real sampleRate() const { return _sampleRate; }
  void setSampleRate(Real sampleRate, bool callSetup = true ) { _sampleRate = sampleRate; if ( callSetup ) setup(); }

  int memorySize() const { return _memorySize; }
  void setMemorySize(int memorySize, bool callSetup = true ) { _memorySize = memorySize; if ( callSetup ) setup(); }

  /**
     Returns the size of the FFT to be performed.  The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const { return _fftSize; };

  /**
     Specifies the @a size of the FFT to be performed.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 will perform faster.
     
     @sa fftSize()
  */
  void setFftSize( int size, bool callSetup = true ) { _fftSize = size; if ( callSetup ) setup(); }

};

#endif  // VOICEACTIVITYDETECTION_H
