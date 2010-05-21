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

#include "BarkBands.h"

#include <vector>

using namespace std;
using namespace Eigen;

BarkBands::BarkBands(int lowBand, int highBand, Real sampleRate, int fftSize) 
{
  
  LOUDIA_DEBUG("BARKBANDS: Constructor lowBand: " << _lowBand << 
        ", highBand: " << _highBand << 
        ", sampleRate: " << _sampleRate << 
        ", fftSize: " << _fftSize);

  if ( lowBand >= highBand ) {
    // Throw an exception, highBand must be higher than lowBand
  }
  
  setLowBand( lowBand, false );
  setHighBand( highBand, false );
  setSampleRate( sampleRate, false );
  setFftSize( fftSize, false );
  
  setup();
  
  LOUDIA_DEBUG("BARKBANDS: Constructed");
  
}

void BarkBands::setup(){
  LOUDIA_DEBUG("BARKBANDS: Setting up...");

  // In some cases the first boundary is set to 0
  MatrixXR startFreqs(25, 1);
  startFreqs << 20, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500;

  MatrixXR centerFreqs(24, 1);
  centerFreqs << 60, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500;
  
  MatrixXI startBins = ((startFreqs.block(_lowBand, 0, _highBand - _lowBand + 2, 1) / _sampleRate) * _fftSize).cast<int>();
  
  std::vector<MatrixXR> weights;
  for (int i = 0; i < startBins.rows() - 1; i++) {
    MatrixXR bandWeights = MatrixXR::Ones(startBins(i+1) - startBins(i), 1);
    weights.push_back(bandWeights);
  }

  _bands.setStartsWeights(startBins.block(0, 0, _highBand - _lowBand + 1, 1), weights);

  LOUDIA_DEBUG("BARKBANDS: Finished set up...");
}

void BarkBands::process(const MatrixXR& spectrum, MatrixXR* bands) {  
  _bands.process(spectrum, bands);
}

void BarkBands::reset(){
  // Initial values
  _bands.reset();
}

void BarkBands::starts(MatrixXI* result) const {
  return _bands.starts( result );
}

void BarkBands::centers(MatrixXR* result) const {
  (*result) = _centersLinear;
}

std::vector<MatrixXR> BarkBands::weights() const {
  return _bands.weights();
}

void BarkBands::bandWeights(int band, MatrixXR* bandWeights) const {
  return _bands.bandWeights( band, bandWeights );
}

int BarkBands::bands() const {
  return _bands.bands();
}

Real BarkBands::lowBand() const{
  return _lowBand;
}
  
void BarkBands::setLowBand( Real band, bool callSetup ){
  _lowBand = band;
  if ( callSetup ) setup();
}

Real BarkBands::highBand() const{
  return _highBand;
}
  
void BarkBands::setHighBand( Real band, bool callSetup ){
  _highBand = band;
  if ( callSetup ) setup();
}

Real BarkBands::sampleRate() const{
  return _sampleRate;
}
  
void BarkBands::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

int BarkBands::fftSize() const{
  return _fftSize;
}

void BarkBands::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}
