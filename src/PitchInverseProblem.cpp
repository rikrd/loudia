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

#include "PitchInverseProblem.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

PitchInverseProblem::PitchInverseProblem(int fftSize, Real lowFrequency, Real highFrequency, Real samplerate, int pitchCount, int harmonicCount, int frequencyCandidateCount, Real peakWidth)
{
  DEBUG("PITCHINVERSEPROBLEM: Construction fftSize: " << fftSize
        << " samplerate: " << samplerate
        << " lowFrequency: " << lowFrequency
        << " highFrequency: " << highFrequency
        << " pitchCount: " << pitchCount
        << " harmonicCount: " << harmonicCount
        << " frequencyCandidateCount: " << frequencyCandidateCount
        << " peakWidth: " << peakWidth);


  setFftSize( fftSize, false );
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setPitchCount( pitchCount, false );
  setHarmonicCount( harmonicCount, false );
  setFrequencyCandidateCount( frequencyCandidateCount, false );
  setPeakWidth( peakWidth, false );
  setSamplerate( samplerate, false );
  setup();
}

PitchInverseProblem::~PitchInverseProblem(){
}

void PitchInverseProblem::setup(){
  DEBUG("PITCHINVERSEPROBLEM: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;


  _peak.setPeakCount( _pitchCount, false );
  _peak.setSortMethod( PeakDetection::BYMAGNITUDE, false );
  _peak.setup();

  _peakInterp.setup();

  
  int frequencyCount = -1 ? _halfSize : _frequencyCandidateCount;

  _regularisation = 2.0;

  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  _inharmonicity = 0.0;

  MatrixXR freqs;
  range(_lowFrequency, _highFrequency, frequencyCount, &freqs);

  _projectionMatrix.resize(_halfSize, freqs.cols());  // We add one that will be the noise component
  _projectionMatrix.setZero();

  DEBUG("PITCHINVERSEPROBLEM: Setting up the projection matrix...");
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int col = 0; col < _projectionMatrix.cols(); col++ ) {
      Real f = freqs(0, col);

      for ( int harmonicIndex = 1; harmonicIndex < _harmonicCount + 1; harmonicIndex++ ) {
        Real mu = harmonicPosition(f, _lowFrequency, _highFrequency, harmonicIndex);
        Real a = harmonicWeight(f, _lowFrequency, _highFrequency, harmonicIndex);
        Real fi = harmonicSpread(f, _lowFrequency, _highFrequency, harmonicIndex);

        _projectionMatrix(row, col) += a * gaussian(row, mu, fi);
      }
    }
  }
  
  /*
  MatrixXR gauss;
  MatrixXR mu;
  MatrixXR a;
  MatrixXR fi;

  //MatrixXR x;
  //range(0, _frequencyCandidateCount, frequencyCount, _projectionMatrix.rows(), &x);
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int harmonicIndex = 1; harmonicIndex < _harmonicCount + 1; harmonicIndex++ ) {
      harmonicPosition(freqs, _lowFrequency, _highFrequency, harmonicIndex, &mu);
      harmonicWeight(freqs, _lowFrequency, _highFrequency, harmonicIndex, &a);
      harmonicSpread(freqs, _lowFrequency, _highFrequency, harmonicIndex, &fi);

      gaussian(row, mu, fi, &gauss);
      
      _projectionMatrix.row(row) += a.cwise() * gauss;
    }
  }
  */

  MatrixXR sourceWeight = MatrixXR::Identity( frequencyCount, frequencyCount );
  MatrixXR targetWeight = MatrixXR::Identity( _halfSize, _halfSize );

  MatrixXR invSourceWeight = ( sourceWeight ).inverse();

  DEBUG("PITCHINVERSEPROBLEM: Setting up the inversion...");
  // A = W^{-1} K^t [ K W^{-1} K^t + \lambda * I_N ]^{+}
  _inverseProjectionMatrix = invSourceWeight * _projectionMatrix.transpose() * ( (_projectionMatrix * invSourceWeight * _projectionMatrix.transpose()) + (_regularisation * MatrixXR::Identity( _halfSize, _halfSize )) ).inverse();

  DEBUG("PITCHINVERSEPROBLEM: Setting up the peak detector and interpolator...");
  _peak.setup();
  _peakInterp.setup();

  reset();

  DEBUG("PITCHINVERSEPROBLEM: Finished setup.");
}


void PitchInverseProblem::harmonicWeight(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result){
  (*result) = MatrixXR::Constant(f.rows(), f.cols(), (fMax + _alpha) / ((harmonicIndex * fMin) + _beta));
}

void PitchInverseProblem::harmonicPosition(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result){
  (*result) = (harmonicIndex * f * sqrt(1.0 + (pow(harmonicIndex, 2.0) - 1.0) * _inharmonicity)) * (Real)_fftSize / (2.0 * (Real)_samplerate);
}

void PitchInverseProblem::harmonicSpread(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result){
  (*result) = MatrixXR::Constant(f.rows(), f.cols(), _peakWidth);
}

Real PitchInverseProblem::harmonicWeight(Real f, Real fMin, Real fMax, int harmonicIndex){
  //return ((_samplerate / tLow) + _alpha) / ((harmonicIndex * _samplerate / tUp) + _beta);
  return ((harmonicIndex * f) + _beta) / ((harmonicIndex * f) + _alpha);
  //return 1.0;
}

Real PitchInverseProblem::harmonicPosition(Real f, Real fMin, Real fMax, int harmonicIndex){
  return (harmonicIndex * f * sqrt(1.0 + (pow(harmonicIndex, 2.0) - 1.0) * _inharmonicity)) * (Real)_fftSize / (2.0 * (Real)_samplerate);
}

Real PitchInverseProblem::harmonicSpread(Real f, Real fMin, Real fMax, int harmonicIndex){
  // TODO: change this by a spread function which might or might not change with the position
  //       or other things such as the chirp rate or inharmonicity error
  return _peakWidth;
}


void PitchInverseProblem::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* freqs){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, _pitchCount );
  (*saliencies).resize( rows, _pitchCount );
  (*freqs).resize( rows, _projectionMatrix.cols() );

  (*pitches).setZero();
  (*saliencies).setZero();

  for ( int row = 0; row < rows; row++ ) {
    DEBUG("PITCHINVERSEPROBLEM: Matrix multiplication");
    DEBUG("PITCHINVERSEPROBLEM: _inverseProjectionMatrix: " << _inverseProjectionMatrix.rows() << ", " << _inverseProjectionMatrix.cols() );
    DEBUG("PITCHINVERSEPROBLEM: spectrum: " << spectrum.rows() << ", " << spectrum.cols() );
    (*freqs).row( row ) = _inverseProjectionMatrix * spectrum.row( row ).transpose();
  }

  DEBUG("PITCHINVERSEPROBLEM: Find peaks");
  
  _peak.process((*freqs),
                pitches, saliencies);

  DEBUG("PITCHINVERSEPROBLEM: Interpolate peaks");
  
  _peakInterp.process((*freqs), (*pitches), (*saliencies),
                      pitches, saliencies);

  DEBUG("PITCHINVERSEPROBLEM: Setting the pitches");
  
  (*pitches) = (((_highFrequency - _lowFrequency) / _frequencyCandidateCount) * (*pitches)).cwise() + _lowFrequency;
  
}

void PitchInverseProblem::reset(){
  // Initial values
  DEBUG("PITCHINVERSEPROBLEM: Resetting...");
  _peak.reset();
  _peakInterp.reset();
  
}

void PitchInverseProblem::projectionMatrix(MatrixXR* matrix) const {
  (*matrix) = _projectionMatrix;
}

Real PitchInverseProblem::lowFrequency() const{
  return _lowFrequency;
}
  
void PitchInverseProblem::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real PitchInverseProblem::highFrequency() const{
  return _highFrequency;
}
  
void PitchInverseProblem::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real PitchInverseProblem::samplerate() const{
  return _samplerate;
}
  
void PitchInverseProblem::setSamplerate( Real frequency, bool callSetup ){
  _samplerate = frequency;
  if ( callSetup ) setup();
}

int PitchInverseProblem::fftSize() const{
  return _fftSize;
}

void PitchInverseProblem::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

int PitchInverseProblem::peakWidth() const{
  return _peakWidth;
}

void PitchInverseProblem::setPeakWidth( int width, bool callSetup ) {
  _peakWidth = width;
  if ( callSetup ) setup();
}

int PitchInverseProblem::frequencyCandidateCount() const{
  return _frequencyCandidateCount;
}

void PitchInverseProblem::setFrequencyCandidateCount( int count, bool callSetup ) {
  _frequencyCandidateCount = count;
  if ( callSetup ) setup();
}

int PitchInverseProblem::pitchCount() const{
  return _pitchCount;
}

void PitchInverseProblem::setPitchCount( int count, bool callSetup ) {
  _pitchCount = count;
  if ( callSetup ) setup();
}

int PitchInverseProblem::harmonicCount() const{
  return _harmonicCount;
}

void PitchInverseProblem::setHarmonicCount( int count, bool callSetup ) {
  _harmonicCount = count;
  if ( callSetup ) setup();
}
