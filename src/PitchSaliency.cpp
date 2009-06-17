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

#include "PitchSaliency.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

PitchSaliency::PitchSaliency(int fftSize, Real f0, Real f1, Real sampleRate, Real fPrec, int numHarmonics) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _f0( f0 ),
  _f1( f1 ),
  _fPrec( fPrec ),
  _numHarmonics( numHarmonics ),
  _sampleRate( sampleRate )
{
  LOUDIA_DEBUG("PITCHSALIENCY: Construction fftSize: " << _fftSize
        << " sampleRate: " << _sampleRate
        << " f0: " << _f0
        << " f1: " << _f1
        << " fPrec: " << _fPrec
        << " numHarmonics: " << _numHarmonics );

  setup();
}

PitchSaliency::~PitchSaliency(){}

void PitchSaliency::setup(){
  LOUDIA_DEBUG("PITCHSALIENCY: Setting up...");
  
  _halfSize = ( _fftSize / 2 ) + 1;

  _tMax = _sampleRate / _f0;
  _tMin = _sampleRate / _f1;
  
  _tPrec = _fPrec;


  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  
  reset();

  LOUDIA_DEBUG("PITCHSALIENCY: Finished setup.");
}

Real PitchSaliency::harmonicWeight(Real period, Real tLow, Real tUp, int harmonicIndex){
  return ((_sampleRate / tLow) + _alpha) / ((harmonicIndex * _sampleRate / tUp) + _beta);
}

Real PitchSaliency::saliency(Real period, Real deltaPeriod, Real tLow, Real tUp, const MatrixXR& spectrum){
  const int cols = spectrum.cols();
  Real sum = 0.0;
  
  for ( int m = 1; m < _numHarmonics; m++ ) {
    
    int begin = (int)round(m * _fftSize / (period + (deltaPeriod / 2.0)));
    int end = min((int)round(m * _fftSize / (period - (deltaPeriod / 2.0))), cols - 1);

    if (begin < end) sum += harmonicWeight(period, tLow, tUp, m) * spectrum.block(0, begin, 1, end - begin).maxCoeff();
  }

  return sum;
}

void PitchSaliency::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, 1 );
  (*saliencies).resize( rows, 1 );
  
  for ( int row = 0; row < rows; row++ ) {

    Real tLow = _tMin;
    Real tUp = _tMax;
    Real sal;

    Real tLowBest = tLow;
    Real tUpBest = tUp;
    Real salBest;
  
    Real period;
    Real deltaPeriod;

    while ( ( tUp - tLow ) > _tPrec ) {
      // Split the best block and compute new limits
      tLow = (tLowBest + tUpBest) / 2.0;
      tUp = tUpBest;
      tUpBest = tLow;
      
      // Compute new saliences for the new blocks
      period = (tLowBest + tUpBest) / 2.0;
      deltaPeriod = tUpBest - tLowBest;
      salBest = saliency(period, deltaPeriod, tLowBest, tUpBest, spectrum.row( row ));

      period = (tLow + tUp) / 2.0;
      deltaPeriod = tUp - tLow;
      sal = saliency(period, deltaPeriod, tLow, tUp, spectrum.row( row ));

      if (sal > salBest) {
        tLowBest = tLow;
        tUpBest = tUp;
        salBest = sal;
      }
    }

    period = (tLowBest + tUpBest) / 2.0;
    deltaPeriod = tUpBest - tLowBest;

    (*pitches)(row, 0) = _sampleRate / period;
    (*saliencies)(row, 0) = saliency(period, deltaPeriod, tLowBest, tUpBest, spectrum.row( row ));
  }
}

void PitchSaliency::reset(){
  // Initial values

}
