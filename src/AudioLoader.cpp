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
#include <libavutil/mathematics.h>

using namespace std;

#define MAX_AUDIO_FRAME_SIZE 192000

AudioLoader::AudioLoader( const std::string& filename, const int frameSize, int channel, Real loadDuration ) :
  _inputFrameSize(MAX_AUDIO_FRAME_SIZE),
  _formatContext(0),
  _audioCodecContext(0),
  _audioStream(-1),
  _sampleRate(0),
  _channelCount(0),
  _buffer(0),
  _audioBuffer(0),
  _audioPacketData(0),
  _finished(false),
  _sizeRead(0)
{
  LOUDIA_DEBUG("AUDIOLOADER: Constructing...");

  av_register_all();

  LOUDIA_DEBUG("AUDIOLOADER: Registered all...");

  setFilename(filename, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the filename...");

  setFrameSize(frameSize, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the frame size...");

  setChannel(channel, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the channel...");

  setLoadDuration(loadDuration, false);

  LOUDIA_DEBUG("AUDIOLOADER: Set the load size...");

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

  // Set the load duration in timebase
  if ( _loadDuration < 0 ) {
    _loadDurationInTimeBase = -1;
  } else {
    //AVStream *stream = _formatContext->streams[_audioStream];
    //_loadDurationInTimeBase = av_rescale(_loadDuration, stream->time_base.den, stream->time_base.num);
  }


  // Create buffer for output audio data
  delete [] _buffer;
  _buffer = new sample_type[_frameSize * _channelCount];

  // Create buffer for input audio data (decoded)
  //_bufferSize = (_inputFrameSize * _channelCount * 3) / 2;
  _bufferSize = _inputFrameSize;
  delete [] _audioBuffer;
  _audioBuffer = (sample_type*)av_malloc(_bufferSize + FF_INPUT_BUFFER_PADDING_SIZE);
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

  case MONOMIX:
    audio->resize(_frameSize, 1);
    audio->setZero();
    for (int i=0, j=0; i < _frameSize; i++, j+=_channelCount) {
      for (int k=0; k < _channelCount; k++){
        (*audio)(i, 0) += scale(_buffer[j+k]) / (Real)_channelCount;
      }
    }
    break;

  default:
    const int channel = _channel%_channelCount;
    audio->resize(_frameSize, 1);
    for (int i=0, j=0; i < _frameSize; i++, j+=_channelCount) {
      (*audio)(i, 0) = scale(_buffer[j+channel]);
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
  }
}

void AudioLoader::loadFile(){
  // Open file
  if (avformat_open_input(&_formatContext, _filename.c_str(), NULL, NULL)!=0) {
    LOUDIA_ERROR("AUDIOLOADER: Could not open file \"" << _filename << "\".  Please set an existing filename using setFilename().");
    return; // Couldn't open file
  }

  // Retrieve stream information
  if (avformat_find_stream_info(_formatContext, NULL)<0) {
    LOUDIA_ERROR("AUDIOLOADER: Could not find stream information in file!");
    return; // Couldn't find stream information
  }

  // Look for an audio stream
  _audioStream = -1;
  for (int i=0; i < (int)_formatContext->nb_streams; i++) {
#if LIBAVCODEC_VERSION_INT >= AV_VERSION_INT(52, 64, 0)
    if (_formatContext->streams[i]->codec->codec_type == AVMEDIA_TYPE_AUDIO) {
#else
    if (_formatContext->streams[i]->codec->codec_type == CODEC_TYPE_AUDIO) {
#endif
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
  if (avcodec_open2(_audioCodecContext, _aCodec, NULL) < 0 ) {
    LOUDIA_ERROR("AUDIOLOADER: Unable to load codec!");
    return;
  }

  _currentTime = 0;
  _loadedDuration = 0;
}

void AudioLoader::closeFile(){
  if (_audioCodecContext) {
    // Close the codec
    avcodec_close(_audioCodecContext);
    _audioCodecContext = 0;
  }

  if (_formatContext) {
    // Close the video file
    avformat_close_input(&_formatContext);
    _formatContext = 0;
  }
}

void AudioLoader::seek( Real time ) {
  if (_audioCodecContext == 0) {
    LOUDIA_ERROR("The algorithm must be set up before hand calling setup().");
  }

  //AVStream *stream = _formatContext->streams[_audioStream];
  //int64_t pts = av_rescale((int64_t)time, stream->time_base.den, stream->time_base.num);
  //if ( av_seek_frame(_formatContext, _audioStream, pts, AVSEEK_FLAG_ANY & AVSEEK_FLAG_BACKWARD) < 0 )
  //  LOUDIA_WARNING("AUDIOLOADER: Unable to seek.");

  _loadedDuration = 0;
}

Real AudioLoader::fileProgress() const {
  // TODO: check if there is a more correct way
    
  Real fileSize = _formatContext->pb ? avio_size(_formatContext->pb) : -1;
  
  if (fileSize == 0) return -1;

  return (Real)_sizeRead / fileSize;
}

Real AudioLoader::loadProgress() const {
  if (_audioCodecContext == 0) return 0;

  AVStream *stream = _formatContext->streams[_audioStream];

  Real totalDuration;

  Real currentDuration = (Real)_loadedDuration * stream->time_base.num / stream->time_base.den;

  if (_loadDurationInTimeBase < 0) {
    totalDuration = (Real)stream->duration * stream->time_base.num / stream->time_base.den;
  }else{
    totalDuration = (Real)_loadDurationInTimeBase * stream->time_base.num / stream->time_base.den;
  }

  return currentDuration / totalDuration;

}

Real AudioLoader::currentTime() const {
  if (_audioCodecContext == 0) return 0;

  AVStream *stream = _formatContext->streams[_audioStream];

  return _currentTime * stream->time_base.num / stream->time_base.den;
}

Real AudioLoader::totalTime() const {
  if (_audioCodecContext == 0) return 0;

  AVStream *stream = _formatContext->streams[_audioStream];

  return (Real)stream->duration * (Real)stream->time_base.num / (Real)stream->time_base.den;
}

bool AudioLoader::nextPacket(){
  while(av_read_frame(_formatContext, &_packet) >= 0) {
    // Is this a packet from the audio stream?
    if( _packet.stream_index == _audioStream && ((_loadDurationInTimeBase < 0) || (_loadedDuration < _loadDurationInTimeBase)) ) {
      _sizeRead = _packet.pos;
      _currentTime = _packet.pts;
      _loadedDuration += _packet.duration;
      return true;

    } else {
      // This packet does not correspond to the stream
      av_free_packet(&_packet);
    }

  }

  // There were no packets left
  return false;
}

int avcodec_decode_audio3(AVCodecContext *avctx, int16_t *samples,
                          int *frame_size_ptr,
                          AVPacket *avpkt)
{
    AVFrame frame;
    int ret, got_frame = 0;

    if (avctx->get_buffer != avcodec_default_get_buffer) {
        av_log(avctx, AV_LOG_ERROR, "Custom get_buffer() for use with"
               "avcodec_decode_audio3() detected. Overriding with avcodec_default_get_buffer\n");
        av_log(avctx, AV_LOG_ERROR, "Please port your application to "
               "avcodec_decode_audio4()\n");
        avctx->get_buffer = avcodec_default_get_buffer;
        avctx->release_buffer = avcodec_default_release_buffer;
    }

    ret = avcodec_decode_audio4(avctx, &frame, &got_frame, avpkt);

    if (ret >= 0 && got_frame) {
        int ch, plane_size;
        int planar = av_sample_fmt_is_planar(avctx->sample_fmt);
        int data_size = av_samples_get_buffer_size(&plane_size, avctx->channels,
                                                   frame.nb_samples,
                                                   avctx->sample_fmt, 1);
        if (*frame_size_ptr < data_size) {
            av_log(avctx, AV_LOG_ERROR, "output buffer size is too small for "
                   "the current frame (%d < %d)\n", *frame_size_ptr, data_size);
            return AVERROR(EINVAL);
        }

        memcpy(samples, frame.extended_data[0], plane_size);

        if (planar && avctx->channels > 1) {
            uint8_t *out = ((uint8_t *)samples) + plane_size;
            for (ch = 1; ch < avctx->channels; ch++) {
                memcpy(out, frame.extended_data[ch], plane_size);
                out += plane_size;
            }
        }
        *frame_size_ptr = data_size;
    } else {
        *frame_size_ptr = 0;
    }
    return ret;
}

int AudioLoader::decodePacket(sample_type* _audioBuffer, int _bufferSize){
  if (_audioCodecContext == 0) {
    LOUDIA_ERROR("The algorithm must be set up before hand calling setup().");
  }

  int decodedSize, dataSize;

  for(;;) {
    while(_audioPacketSize > 0) {
      dataSize = _bufferSize;
#if LIBAVCODEC_VERSION_INT >= AV_VERSION_INT(52, 64, 0)
      AVPacket avpkt;
      av_init_packet(&avpkt);
      avpkt.data = _audioPacketData;
      avpkt.size = _audioPacketSize;
      decodedSize = avcodec_decode_audio3(_audioCodecContext,
                                   _audioBuffer, &dataSize,
                                   &avpkt);
#else

      decodedSize = avcodec_decode_audio2(_audioCodecContext,
                                   _audioBuffer, &dataSize,
                                   _audioPacketData, _audioPacketSize);


#endif

      if(decodedSize < 0) {
      /* if error, skip frame */
        LOUDIA_WARNING("FFMPEG decoding error. Skipping frame.");
        _audioPacketSize = 0;
        break;
      }

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


