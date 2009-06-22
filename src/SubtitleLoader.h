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

#ifndef SUBTITLELOADER_H
#define SUBTITLELOADER_H

#include "Typedefs.h"
#include "Debug.h"

#include <string> 
#include <vector>
#include <cmath> 

#include <boost/regex.hpp>


class SubtitleLoader {
private:
  std::vector<Real*> _timestamps;
  std::string _filename;
  
  Real _offset;
  
  Real _sampleRate;
  Real _loadDuration;
  
  int _streamIndex;
  int _subtitleIndex;
  
  int _frameSize;
  int _hopSize;
  
  bool _finished;

  boost::regex _expression;

  Real seconds(int hours, int mins, int secs, int msecs);

  Real* parseTime(const std::string& str);
  
public:
  //Subtitle loader class
  SubtitleLoader(const std::string& filename = "", const Real sampleRate = 44100.0, const int frameSize = 1024, int hopSize = -1, Real loadTime = -1.0);
  ~SubtitleLoader();

  // Prepare a frame for output (a frame can be filled by a fraction, one or several packets)
  void process( MatrixXR* subtitle );

  void setup();

  void setFilename(const char* filename, const bool callSetup = true) { _filename = std::string(filename); if ( callSetup ) setup(); };
  void setFilename(const std::string& filename, const bool callSetup = true) { _filename = filename; if ( callSetup ) setup(); };
  const std::string& filename() const { return _filename; };

  void setFrameSize( int size, const bool callSetup = true ) { _frameSize = size; if ( callSetup ) setup(); };
  int frameSize() const { return _frameSize; };

  void setLoadDuration( Real duration, const bool callSetup = true ) { _loadDuration = duration; if ( callSetup ) setup(); };
  Real loadDuration() const { return _loadDuration; };

  /**
     Returns the hop size of the algorithm.
     
     By default it is 1024.
  */
  int hopSize() const {
    if ( _hopSize <= 0 ) {
      return std::max(_frameSize / 2, 1);
    }
    
    return _hopSize;
  };
  
  /**
     Specifies the hop @a size of the algorithm.
  */
  void setHopSize( int size, bool callSetup = true ) { _hopSize = size; if (callSetup) setup(); };

  /**
     Return the sampleRate frequency of the input signal.
     The default is 44100.0.

     @sa setSampleRate
  */  
  Real sampleRate() const { return _sampleRate; };  

  /**
     Specifies the sampleRate @a frequency of the input signal.
     
     @sa sampleRate
  */
  void setSampleRate( Real frequency, bool callSetup = true ) {_sampleRate = frequency; if (callSetup) setup();}

  void seek( Real time );
  
  bool finished() const { return _finished; };
  bool isFinished() const { return _finished; };
};

#endif //#ifndef SUBTITLELOADER_H
