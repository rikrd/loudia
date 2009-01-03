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

#include <Eigen/Core>
#include <Eigen/Array>

#include <iostream>
#include <cmath>
#include <vector>

#include "chebyshev.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

Chebyshev::Chebyshev(int channels, int order, Real rippleDB, Real samplerate) : _filter(channels)
{
  DEBUG("CHEBYSHEV: Constructor order: " << order << ", rippleDB: " << rippleDB << ", samplerate: " << samplerate);

  if ( order < 1 ) {
    // Throw an exception
  }

  if ( samplerate <= 0 ) {
    // Throw an exception
  }
  
  _channels = channels;
  _order = order;
  _rippleDB = rippleDB;
  _samplerate = samplerate;

  setup();
  DEBUG("CHEBYSHEV: Constructed");
}

void Chebyshev::setup(){
  DEBUG("CHEBYSHEV: Setting up...");
  
  MatrixXC zeros(1,1);
  zeros.setZero();
  
  Real eps = sqrt(pow(10, (0.1 * _rippleDB) - 1.0));

  MatrixXR n(1, _order + 1);
  for ( int i = 0; i < n.cols(); i++ ) {
    n(0, i) = i;
  }
  
  Real mu = 1.0 / _order * log((1.0 + sqrt( 1 + eps * eps)) / eps);

  MatrixXC theta = ((n * 2).cwise() - 1.0) / _order * M_PI / 2.0;

  MatrixXC poles = -sinh(mu) * theta.cwise().sin() + Complex(0, 1) * cosh(mu) * theta.cwise().cos();

  Complex gainComplex = 1.0;
  for ( int i = 0; i < poles.cols(); i++ ) {
    gainComplex *= -poles(0, i);
  }

  Real gain = gainComplex.real();
  
  if ( _order % 2 == 0 ) {
    gain = gain / sqrt((1 + eps * eps));
  }
  
  DEBUG("CHEBYSHEV: zeros:" << zeros );
  DEBUG("CHEBYSHEV: poles:" << poles );
  DEBUG("CHEBYSHEV: gain:" << gain );

  DEBUG("CHEBYSHEV: Finished set up...");
}

void Chebyshev::process(MatrixXR samples, MatrixXR* filtered) {
  _filter.process(samples, filtered);
}

void Chebyshev::reset(){
  // Initial values
  _filter.reset();
}
