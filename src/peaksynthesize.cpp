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

#include "typedefs.h"
#include "debug.h"

#include "peaksynthesize.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

PeakSynthesize::PeakSynthesize(int fftSize) {
  DEBUG("PEAKSYNTHESIZE: Constructor");
  
  _fftSize = fftSize;

  setup();
  
  DEBUG("PEAKSYNTHESIZE: Constructed");
}

PeakSynthesize::~PeakSynthesize() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakSynthesize::setup(){
  // Prepare the buffers
  DEBUG("PEAKSYNTHESIZE: Setting up...");

  reset();

  DEBUG("PEAKSYNTHESIZE: Finished set up...");
}


void PeakSynthesize::process(const MatrixXR& trajPositions, const MatrixXR& trajMagnitudes,
                             MatrixXR* spectrum){
  
  DEBUG("PEAKSYNTHESIZE: Processing");
  
  spectrum->resize(trajPositions.rows(), _fftSize);
  
  for ( int row = 0 ; row < spectrum->rows(); row++ ) {
  
    for ( int i = 0; i < trajPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( trajPositions(row, i) != -1 ){
        
        
      }
    }
  }
  
  DEBUG("PEAKSYNTHESIZE: Finished Processing");
}

void PeakSynthesize::reset(){
  // Initial values
}
