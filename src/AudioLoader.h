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

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

extern "C" {
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#ifdef LOUDIA_OLD_FFMPEG
#include <ffmpeg/avcodec.h>
#include <ffmpeg/avformat.h>
#else
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#endif
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
    ALL                   = -1 /**< All the cahnnels available */,
    MONOMIX               = -2 /**< Mono mix of all the channels */,
    LEFT                  =  0 /**< First channel in the stream (left) */,
    RIGHT                 =  1 /**< Second channel in the stream (right) */
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

  int _bufferSize;
  sample_type* _audioBuffer;
  unsigned int _audioBufferSize;
  unsigned int _audioBufferIndex;

  uint8_t *_audioPacketData;
  int _audioPacketSize;

  bool _finished;

  int _channel;

  int64_t _sizeRead;
  int64_t _currentTime;

  Real _loadDuration;
  int64_t _loadDurationInTimeBase;
  int64_t _loadedDuration;

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
  AudioLoader(const std::string& filename = "", const int frameSize = 1024, int channel = ALL, Real loadTime = -1.0);
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

  void setLoadDuration( Real duration, const bool callSetup = true ) { _loadDuration = duration; if ( callSetup ) setup(); };
  Real loadDuration() const { return _loadDuration; };

  Real loadProgress() const;
  Real fileProgress() const;
  Real currentTime() const;
  Real totalTime() const;

  void seek( Real time );

  bool finished() const { return _finished; };
};

#endif //#ifndef AUDIOLOADER_H
