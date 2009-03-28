/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or   
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

PitchInverseProblem::PitchInverseProblem(int fftSize, Real lowFrequency, Real highFrequency, Real samplerate, int numMaxPitches, int numHarmonics, int numFreqCandidates, Real peakBandwidth) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _lowFrequency( lowFrequency ),
  _highFrequency( highFrequency ),
  _numMaxPitches( numMaxPitches ),
  _numHarmonics( numHarmonics ),
  _numFreqCandidates( numFreqCandidates == -1 ? _halfSize : numFreqCandidates ),
  _peakBandwidth( peakBandwidth ),
  _samplerate( samplerate ),
  _peak(_numMaxPitches, PeakDetection::BYMAGNITUDE ),
  _peakInterp()

{
  DEBUG("PITCHINVERSEPROBLEM: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate
        << " lowFrequency: " << _lowFrequency
        << " highFrequency: " << _highFrequency
        << " numMaxPitches: " << _numMaxPitches
        << " numHarmonics: " << _numHarmonics
        << " numFreqCandidates: " << _numFreqCandidates
        << " peakBandwidth: " << _peakBandwidth);

  setup();
}

PitchInverseProblem::~PitchInverseProblem(){
}

void PitchInverseProblem::setup(){
  DEBUG("PITCHINVERSEPROBLEM: Setting up...");
  
  _regularisation = 2.0;

  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  _inharmonicity = 0.0;

  MatrixXR freqs;
  range(_lowFrequency, _highFrequency, _numFreqCandidates, &freqs);

  _projectionMatrix.resize(_halfSize, freqs.cols());  // We add one that will be the noise component
  _projectionMatrix.setZero();

  DEBUG("PITCHINVERSEPROBLEM: Setting up the projection matrix...");
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int col = 0; col < _projectionMatrix.cols(); col++ ) {
      Real f = freqs(0, col);

      for ( int harmonicIndex = 1; harmonicIndex < _numHarmonics + 1; harmonicIndex++ ) {
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
  //range(0, _numFreqCandidates, _numFreqCandidates, _projectionMatrix.rows(), &x);
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int harmonicIndex = 1; harmonicIndex < _numHarmonics + 1; harmonicIndex++ ) {
      harmonicPosition(freqs, _lowFrequency, _highFrequency, harmonicIndex, &mu);
      harmonicWeight(freqs, _lowFrequency, _highFrequency, harmonicIndex, &a);
      harmonicSpread(freqs, _lowFrequency, _highFrequency, harmonicIndex, &fi);

      gaussian(row, mu, fi, &gauss);
      
      _projectionMatrix.row(row) += a.cwise() * gauss;
    }
  }
  */

  MatrixXR sourceWeight = MatrixXR::Identity( _numFreqCandidates, _numFreqCandidates );
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
  (*result) = MatrixXR::Constant(f.rows(), f.cols(), _peakBandwidth);
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
  return _peakBandwidth;
}


void PitchInverseProblem::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* freqs){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, _numMaxPitches );
  (*saliencies).resize( rows, _numMaxPitches );
  (*freqs).resize( rows, _numFreqCandidates );

  (*pitches).setZero();
  (*saliencies).setZero();

  for ( int row = 0; row < rows; row++ ) {
    DEBUG("PITCHINVERSEPROBLEM: Setting the result");
    (*freqs).row( row ) = _inverseProjectionMatrix * spectrum.row( row ).transpose();
  }

  
  _peak.process((*freqs),
                pitches, saliencies);
  
  _peakInterp.process((*freqs), (*pitches), (*saliencies),
                      pitches, saliencies);

  (*pitches) = (((_highFrequency - _lowFrequency) / _numFreqCandidates) * (*pitches)).cwise() + _lowFrequency;
  
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
