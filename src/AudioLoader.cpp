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

#include "AudioLoader.h"

#include <iostream>
using namespace std;

AudioLoader::AudioLoader( const std::string& filename, const int frameSize, int channel ) :
  _inputFrameSize(AVCODEC_MAX_AUDIO_FRAME_SIZE),
  _formatContext(0),
  _audioCodecContext(0),
  _audioStream(-1),
  _sampleRate(0),
  _channelCount(0),
  _buffer(0),
  _packetNumber(0),
  _frameNumber(0),
  _sampleNumber(0),
  _audioBuffer(0),
  _audioPacketData(0),
  _finished(false)
{
  LOUDIA_DEBUG("AUDIOLOADER: Constructing...");

  av_register_all();

  LOUDIA_DEBUG("AUDIOLOADER: Registered all...");

  setFilename(filename, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the filename...");

  setFrameSize(frameSize, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the framesize...");

  setChannel(channel, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the channel...");

  // If the filename has been passed then we setup
  if ( _filename != "" ){
    setup();
  }
}

AudioLoader::~AudioLoader(){
  closeFile();
  delete [] _buffer;
  delete [] _audioBuffer;
}

void AudioLoader::setup(){
  closeFile();
  loadFile();

  // Create buffer for output audio data
  delete [] _buffer;
  _buffer = new sample_type[_frameSize * _channelCount];

  // Create buffer for input audio data (decoded)
  //_bufferSize = (_inputFrameSize * _channelCount * 3) / 2;
  _bufferSize = _inputFrameSize;// + FF_INPUT_BUFFER_PADDINGSIZE;
  delete [] _audioBuffer;
  _audioBuffer = new sample_type[_bufferSize + FF_INPUT_BUFFER_PADDING_SIZE];
  _audioBufferSize = 0;
  _audioBufferIndex = 0;

  // Create buffer for packet data (encoded)
  delete [] _audioPacketData;
  _audioPacketData = 0;
  _audioPacketSize = 0;

  // TODO: check if it is necessary
  av_init_packet(&_packet);

  _finished = false;
}

void AudioLoader::process(MatrixXR *audio){
  process(_buffer);
  
  switch (_channel) {
  case ALL:  
    audio->resize(_frameSize, _channelCount);

    for (int i=0, j=0; i < _frameSize; i++, j+=_channelCount) {
      for (int k=0; k < _channelCount; k++){
        (*audio)(i, k) = scale(_buffer[j+k]);
      }
    }
    
    break;

  case MIX:
    audio->resize(_frameSize, 1);
    for (int i=0, j=0; i < _frameSize; i++, j+=_channelCount) {
      for (int k=0; k < _channelCount; k++){
        (*audio)(i, 0) += scale(_buffer[j+k]) / (Real)_channelCount;
      }
    }    
    break;

  default:
    audio->resize(_frameSize, 1);
    for (int i=0, j=0; i < _frameSize; i++, j+=_channelCount) {
      (*audio)(i, 0) += scale(_buffer[j+_channel]);
    }
    break;
  }
  
}

void AudioLoader::process(sample_type *audioLR){
  int len1;  // temporal output buffer length in bytes
  int audioSize;  // input buffer (decoded audio samples) length in bytes
  int len = _frameSize * _channelCount * sizeof(sample_type);  // target output buffer length in bytes
  sample_type *stream = audioLR;

  while(len > 0) {
    if(_audioBufferIndex >= _audioBufferSize) {
      /* We have already sent all our data; get more */
      audioSize = decodePacket(_audioBuffer,
                                _inputFrameSize * sizeof(sample_type));

      if(audioSize < 0) {
        /* If error, output silence */
        _audioBufferSize = len * sizeof(sample_type);
	memset(_audioBuffer, 0, _audioBufferSize);
      } else {
	_audioBufferSize = audioSize;
      }

      _audioBufferIndex = 0;
    }

    len1 = _audioBufferSize - _audioBufferIndex;

    if(len1 > len)
      len1 = len;

    memcpy(stream, _audioBuffer + _audioBufferIndex/sizeof(sample_type), len1);

    len -= len1;
    stream += len1/sizeof(sample_type);
    _audioBufferIndex += len1;
    _sampleNumber += len1/sizeof(sample_type);
  }
}

void AudioLoader::loadFile(){
  // Open file
  if (av_open_input_file(&_formatContext, _filename.c_str(), NULL, 0, NULL)!=0) {
    LOUDIA_ERROR("AUDIOLOADER: Could not open file \"" << _filename << "\".  Please set an existing filename using setFilename().");
    return; // Couldn't open file
  }

  // Retrieve stream information
  if (av_find_stream_info(_formatContext)<0) {
    LOUDIA_ERROR("AUDIOLOADER: Could not find stream information in file!");
    return; // Couldn't find stream information
  }

  // Look for an audio stream
  _audioStream = -1;
  for (int i=0; i < (int)_formatContext->nb_streams; i++) {
    if (_formatContext->streams[i]->codec->codec_type == CODEC_TYPE_AUDIO) {
      _audioStream = i;
      break;
    }
  }

  // Check if a stream has been found
  if (_audioStream == -1) {
    LOUDIA_ERROR("AUDIOLOADER: No audio stream was found!");
    return;
  }

  // Get the codec's context from the stream
  _audioCodecContext = _formatContext->streams[_audioStream]->codec;

  // Get the info from the stream
  _sampleRate = _audioCodecContext->sample_rate;
  _channelCount = _audioCodecContext->channels;

  // Find the decoder
  _aCodec = avcodec_find_decoder(_audioCodecContext->codec_id);
  if (!_aCodec) {
    LOUDIA_ERROR("AUDIOLOADER: Unsupported codec!");
    return;
  }

  // Load the decoder
  if (avcodec_open(_audioCodecContext, _aCodec) < 0 ) {
    LOUDIA_ERROR("AUDIOLOADER: Unable to load codec!");
    return;
  }
}

void AudioLoader::closeFile(){
  if (_audioCodecContext) {
    // Close the codec
    avcodec_close(_audioCodecContext);
    _audioCodecContext = 0;
  }

  if (_formatContext) {
    // Close the video file
    av_close_input_file(_formatContext);
    _formatContext = 0;
  }
}


bool AudioLoader::nextPacket(){
  while(av_read_frame(_formatContext, &_packet) >= 0) {
    // Is this a packet from the audio stream?
    if( _packet.stream_index == _audioStream ) {
      // We found another packet corresponding to the stream
      _packetNumber++;
      return true;
    } else {
      // This packet does not correspond to the stream
      av_free_packet(&_packet);
    }

  }

  // There were no packets left
  return false;
}

int AudioLoader::decodePacket(sample_type* _audioBuffer, int _bufferSize){
  if (_audioCodecContext == 0) {
    LOUDIA_ERROR("The algorithm must be set up before hand calling setup().");
  }

  int decodedSize, dataSize;

  for(;;) {
    while(_audioPacketSize > 0) {
      dataSize = _bufferSize;

      decodedSize = avcodec_decode_audio2(_audioCodecContext, _audioBuffer, &dataSize,
                                   _audioPacketData, _audioPacketSize);
      if(decodedSize < 0) {
	/* if error, skip frame */
        LOUDIA_WARNING("FFMPEG decoding error. Skipping frame.");
	_audioPacketSize = 0;
	break;
      }

      _frameNumber += decodedSize/sizeof(sample_type);

      _audioPacketData += decodedSize;
      _audioPacketSize -= decodedSize;

      if(dataSize <= 0) {
	/* No data yet, get more frames */
	continue;
      }

      /* We have data, return it and come back for more later */
      return dataSize;
    }

    av_free_packet(&_packet);

    _finished = !nextPacket();
    if (_finished) {
      // We have finished
      return -1;
    }

    _audioPacketData = _packet.data;
    _audioPacketSize = _packet.size;
  }
}


