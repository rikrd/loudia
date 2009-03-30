/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#include "PeakSynthesize.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

PeakSynthesize::PeakSynthesize(int windowSize, int fftSize, Window::WindowType windowType) :
  _windowSize( windowSize ),
  _windowType( windowType ),
  _fftSize( fftSize )

{
  DEBUG("PEAKSYNTHESIZE: Constructor");
  
  setup();
  
  DEBUG("PEAKSYNTHESIZE: Constructed");
}

PeakSynthesize::~PeakSynthesize() {}


void PeakSynthesize::setup(){
  // Prepare the buffers
  DEBUG("PEAKSYNTHESIZE: Setting up...");

  reset();

  DEBUG("PEAKSYNTHESIZE: Finished set up...");
}


void PeakSynthesize::process(const MatrixXR& trajPositions, const MatrixXR& trajMagnitudes,
                             MatrixXR* spectrum){
  
  DEBUG("PEAKSYNTHESIZE: Processing");
  
  spectrum->resize(trajPositions.rows(), (int)ceil(_fftSize/2.0));
  spectrum->setZero();
  
  MatrixXR trajMags;
  dbToMag(trajMagnitudes, &trajMags);
  
  for ( int row = 0 ; row < spectrum->rows(); row++ ) {
  
    for ( int i = 0; i < trajPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( trajPositions(row, i) != -1 ){
        MatrixXR windowTransform;
        int begin, end;

        switch(_windowType){

        case Window::RECTANGULAR:
          // TODO: Implement this window transform

        case Window::HANN:
        case Window::HANNING:
          hannTransform(trajPositions(row, i), trajMags(row, i), 
                        _windowSize, _fftSize, 
                        &windowTransform, &begin, &end);
          break;
          
        case Window::HAMMING:
          hammingTransform(trajPositions(row, i), trajMags(row, i), 
                           _windowSize, _fftSize, 
                           &windowTransform, &begin, &end);
          break;

        default:
          DEBUG("ERROR: Unknown type of window");
          // Throw ValueError unknown window type
          break;

        }
        
        spectrum->block(0, begin, 1, windowTransform.cols()) += windowTransform.row(0);        
      }
    }
  }
  
  DEBUG("PEAKSYNTHESIZE: Finished Processing");
}

void PeakSynthesize::reset(){
  // Initial values
}
