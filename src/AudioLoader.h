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

#ifndef AUDIOLOADER_H
#define AUDIOLOADER_H

#include "Typedefs.h"
#include "Debug.h"
#include <string>

extern "C" {
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
}

typedef int16_t sample_type;

class AudioLoader {
public:
  /**
     @enum Channel
     @brief Specifies the special channels such as ALL and MIX.
     
     @sa channel
  */
  enum Channel {    
    ALL                   = -1 /**< Hold the value until the next sample */,
    MIX                   = -2 /**< Fastest cardinal sine method */
  };
  
private:
  const int _inputFrameSize;
  int _frameSize;
  std::string _filename;

  AVFormatContext *_formatContext;
  AVCodecContext *_audioCodecContext;
  AVCodec *_aCodec;
  AVPacket _packet;
  int _audioStream;

  int _sampleRate;
  int _channelCount;
  
  sample_type *_buffer;

  int _packetNumber;
  int _frameNumber;
  int _sampleNumber;

  int _bufferSize;
  sample_type* _audioBuffer;
  unsigned int _audioBufferSize;
  unsigned int _audioBufferIndex;

  uint8_t *_audioPacketData;
  int _audioPacketSize;

  bool _finished;

  int _channel;

  int _sizeRead;

  void loadFile();
  void closeFile();

  // Get the next packet in the audio stream
  bool nextPacket();

  // Decode the current packet into a buffer
  int decodePacket( sample_type *buffer, int size );

  inline Real scale(sample_type value) {
    return value / (Real)32767;
  }

  void process( sample_type* audio );

public:
  //Audio loader class
  AudioLoader(const std::string& filename = "", const int frameSize = 1024, int channel = ALL);
  ~AudioLoader();

  // Prepare a frame for output (a frame can be filled by a fraction, one or several packets)
  void process( MatrixXR* audio );

  void setup();

  int sampleRate() const { return _sampleRate; };
  int channelCount() const { return _channelCount; };
  bool isFinished() const { return _finished; };

  void setFilename(const char* filename, const bool callSetup = true) { _filename = std::string(filename); if ( callSetup ) setup(); };
  void setFilename(const std::string& filename, const bool callSetup = true) { _filename = filename; if ( callSetup ) setup(); };
  const std::string& filename() const { return _filename; };

  void setFrameSize( int size, const bool callSetup = true ) { _frameSize = size; if ( callSetup ) setup(); };
  int frameSize() const { return _frameSize; };

  void setChannel( int channel, const bool callSetup = true ) { _channel = channel; if ( callSetup ) setup(); };
  int channel() const { return _channel; };  

  float progress() const;
  bool finished() const { return _finished; };
};

#endif //#ifndef AUDIOLOADER_H
