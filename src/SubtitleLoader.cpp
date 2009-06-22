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

#include "SubtitleLoader.h"

#include <iostream>
#include <fstream>
using namespace std;

SubtitleLoader::SubtitleLoader( const std::string& filename, const Real sampleRate, const int frameSize, const int hopSize, Real loadDuration ) :
  _offset(0.0),
  _expression("[0-9][0-9]:[0-9][0-9]:[0-9][0-9],[0-9][0-9][0-9] --> "
              "[0-9][0-9]:[0-9][0-9]:[0-9][0-9],[0-9][0-9][0-9]")
{
  LOUDIA_DEBUG("SUBTITLELOADER: Constructing...");

  setFilename(filename, false);

  LOUDIA_DEBUG("SUBTITLELOADER: Set the filename...");

  setFrameSize(frameSize, false);

  LOUDIA_DEBUG("SUBTITLELOADER: Set the frame size...");

  setHopSize(hopSize, false);

  LOUDIA_DEBUG("SUBTITLELOADER: Set the hop size...");

  setSampleRate(sampleRate, false);

  LOUDIA_DEBUG("SUBTITLELOADER: Set the sample rate...");

  setLoadDuration(loadDuration, false);

  LOUDIA_DEBUG("SUBTITLELOADER: Set the load size...");

  // If the filename has been passed then we setup
  if ( _filename != "" ){
    setup();
  }
}

SubtitleLoader::~SubtitleLoader(){
  for ( int i = 0; i<(int)_timestamps.size(); i++) {
    delete [] _timestamps[i];
  }
  _timestamps.clear();
}

void SubtitleLoader::setup(){
  ifstream subtitleFile;
  subtitleFile.open(_filename.c_str(), ios::in);    // open the streams

  for ( int i = 0; i<(int)_timestamps.size(); i++) {
    delete [] _timestamps[i];
  }
  _timestamps.clear();

  boost::match_results<string::const_iterator> what; 
  boost::match_flag_type flags = boost::match_default;
  string line;

  while(getline(subtitleFile, line)) {
    if (regex_search(line, what, _expression, flags) ) {
      //cout << "--------- " << what[0] << " ---------" << endl;
      Real* times;
      times = parseTime(what[0]);
      _timestamps.push_back(times);
      //cout << "time : " << times[0] << " - " << times[1] << endl;
    }
  }

  subtitleFile.close();
  
  _subtitleIndex = 0;
  _streamIndex = 0;
  _finished = false;
}

void SubtitleLoader::process(MatrixXR *samples){

  samples->resize(_frameSize, 1);
  samples->setZero();
  
  while(true) {    
    // There are no more subtitle
    if (_subtitleIndex >= (int)_timestamps.size()) {
      _finished = true;
      return;
    }
    
    int begin = (_timestamps[_subtitleIndex][0] - _offset) * _sampleRate;
    int end = (_timestamps[_subtitleIndex][1] - _offset) * _sampleRate;
    
    // We are finished since we have passed the load duration
    if ((_loadDuration > 0) && (_streamIndex > (_loadDuration * _sampleRate))) {
      _finished = true;
      return;
    }

    // The current subtitle is after the current frame
    if (begin > (_streamIndex + _frameSize)) break;
    
    // Get the indices of the current subtitle region
    int beginInFrame = max(0, begin - _streamIndex);
    int endInFrame = min(_frameSize, end - _streamIndex);
    
    samples->col(0).segment(beginInFrame, endInFrame - beginInFrame).setOnes();
    
    // The current subtitle is before the current frame
    if (end <= (_streamIndex + _frameSize)) {
      _subtitleIndex++;
    } else {
      break;
    }
  }

  _streamIndex += (_hopSize > 0) ? _hopSize : _frameSize;
  return;
}

void SubtitleLoader::seek( Real time ) {
  _offset = time;
}

Real SubtitleLoader::seconds(int hours, int mins, int secs, int msecs){
  return hours * 3600.0 + mins * 60.0 + secs + msecs/1000.0; 
}

Real* SubtitleLoader::parseTime(const string& str) {
  Real* times = new Real[2];
  int hours1, mins1, secs1, msecs1;
  int hours2, mins2, secs2, msecs2;

  sscanf(str.c_str(), "%02d:%02d:%02d,%03d --> %02d:%02d:%02d,%03d", &hours1, &mins1, &secs1, &msecs1, &hours2, &mins2, &secs2, &msecs2);

  times[0] = seconds(hours1, mins1, secs1, msecs1);
  times[1] = seconds(hours2, mins2, secs2, msecs2);
  
  return times;
}
