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

PitchInverseProblem::PitchInverseProblem(int fftSize, Real f0, Real f1, Real samplerate, Real fPrec, int numHarmonics, int numFreqCandidates, int numMaxPitches, int peakBandwidth) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _f0( f0 ),
  _f1( f1 ),
  _fPrec( fPrec ),
  _numHarmonics( numHarmonics ),
  _numFreqCandidates( numFreqCandidates == -1 ? _halfSize : min( numFreqCandidates, _halfSize ) ),
  _numMaxPitches( numMaxPitches ),
  _peakBandwidth( peakBandwidth ),
  _samplerate( samplerate )
{
  DEBUG("PITCHINVERSEPROBLEM: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate
        << " f0: " << _f0
        << " f1: " << _f1
        << " fPrec: " << _fPrec
        << " numHarmonics: " << _numHarmonics
        << " numFreqCandidates: " << _numFreqCandidates
        << " numMaxPitches: " << _numMaxPitches
        << " peakBandwidth: " << _peakBandwidth);

  setup();
}

PitchInverseProblem::~PitchInverseProblem(){
  delete _inverseProjectionMatrix;
}

void PitchInverseProblem::setup(){
  DEBUG("PITCHINVERSEPROBLEM: Setting up...");

  _tMax = _samplerate / _f0;
  _tMin = _samplerate / _f1;
  
  _tPrec = _fPrec;

  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  _inharmonicity = 0.0;

  MatrixXR freqs;
  range(_f0, _f1, _halfSize, &freqs);

  _projectionMatrix.resize(freqs.cols(), _numFreqCandidates);
  _projectionMatrix.setZero();

  DEBUG("PITCHINVERSEPROBLEM: Setting up the projection matrix...");  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int col = 0; col < _projectionMatrix.cols(); col++ ) {
      for ( int harmonicIndex = 1; harmonicIndex < _numHarmonics+1; harmonicIndex++ ) {
        Real f = freqs(0, row);
        Real mu = harmonicPosition(1.0/f, _tMin, _tMax, harmonicIndex);
        Real a = harmonicWeight(1.0/f, _tMin, _tMax, harmonicIndex);
        Real fi = harmonicSpread(1.0/f, _tMin, _tMax, harmonicIndex);
        
        _projectionMatrix(row, col) += a * gaussian(col, mu, fi);
      }
    }
  }

  DEBUG("PITCHINVERSEPROBLEM: Setting up the LU decomposition...");
  _inverseProjectionMatrix = new LU<MatrixXR>(_projectionMatrix);

  reset();

  DEBUG("PITCHINVERSEPROBLEM: Finished setup.");
}

Real PitchInverseProblem::harmonicWeight(Real period, Real tLow, Real tUp, int harmonicIndex){
  return ((_samplerate / tLow) + _alpha) / ((harmonicIndex * _samplerate / tUp) + _beta);
}

Real PitchInverseProblem::harmonicPosition(Real period, Real tLow, Real tUp, int harmonicIndex){
  return (harmonicIndex / period * sqrt(1.0 + (pow(harmonicIndex, 2.0) - 1.0)*_inharmonicity));
}

Real PitchInverseProblem::harmonicSpread(Real period, Real tLow, Real tUp, int harmonicIndex){
  // TODO: change this by a spread function which might or might not change with the position
  //       or other things such as the chirp rate or inharmonicity error
  return _peakBandwidth;
}


void PitchInverseProblem::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* freqs){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, _numMaxPitches );
  (*saliencies).resize( rows, _numMaxPitches );
  (*freqs).resize( rows, _numFreqCandidates );

  MatrixXR a(1, _numFreqCandidates);

  for ( int row = 0; row < rows; row++ ) {
    DEBUG("PITCHINVERSEPROBLEM: Solving the LU");
    _inverseProjectionMatrix->solve( spectrum.row( row ).transpose(), &a );
    
    DEBUG("PITCHINVERSEPROBLEM: Setting the result");
    (*freqs).row( row ) = a;
  }
}

void PitchInverseProblem::reset(){
  // Initial values

}

void PitchInverseProblem::projectionMatrix(MatrixXR* matrix) const {
  (*matrix) = _projectionMatrix;
}
