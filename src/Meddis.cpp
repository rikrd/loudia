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

#include "Meddis.h"

using namespace std;
using namespace Eigen;

  // Internal fixed constants
const Real Meddis::M = 1.0;
const Real Meddis::A = 5.0;
const Real Meddis::B = 300.0;
const Real Meddis::g = 2000.0;
const Real Meddis::y = 5.05;
const Real Meddis::l = 2500.0;
const Real Meddis::r = 6580.0;
const Real Meddis::x = 66.31;
const Real Meddis::h = 50000.0;

Meddis::Meddis(Real sampleRate, int channels, bool substractSpont) : 
  _sampleRate( sampleRate ),
  _channels( channels ),
  _substractSpont( substractSpont ),
  kt( 1, channels ), 
  spont( 1, channels ), 
  c( 1, channels ), 
  q( 1, channels ), 
  w( 1, channels ) 
{
  LOUDIA_DEBUG("MEDDIS: Constructor sampleRate:" << sampleRate << 
        ", channels:" << channels << 
        ", substractSpont:" << substractSpont);
  
  setup();
}

Meddis::~Meddis() {}

void Meddis::setup(){
  LOUDIA_DEBUG("MEDDIS: Setting up...");

  // Configure the internal constants
  dt = 1./_sampleRate;
  gdt = g*dt;
  ydt = y*dt;
  ldt = l*dt;
  rdt = r*dt;
  xdt = x*dt;

  // Prepare the buffers
  reset();
  LOUDIA_DEBUG("MEDDIS: Set up finished.");
}

void Meddis::process(const MatrixXR& samples, MatrixXR* output){
  // Process will be called with a matrix where columns will be channels
  // and rows will be the time axis
  MatrixXR  row, limitedSt, replenish, eject, loss, reuptake, reprocess;

  (*output).resize(samples.rows(), _channels);

  for (int i = 0; i < samples.rows(); ++i) {
    row = samples.row(i);
  
    limitedSt = (row.cwise() + A).array().clipUnder();

    kt = (limitedSt * gdt).array() / (limitedSt.array() + B);

    replenish = ydt * ((-q).cwise() + M).array().clipUnder();
    eject = kt.cwise() * q;
    loss = ldt * c;
    reuptake = rdt * c;
    reprocess = xdt * w;
      
    q += replenish - eject + reprocess;
    c += eject - loss - reuptake;
    w += reuptake - reprocess;

    // Now iterate through each time slice of the data.  Use the
    // max function to implement the "if (0>" test.
    (*output).row(i) = h * c;
    
    if(_substractSpont){
      (*output).row(i) = ((*output).row(i) - spont.row(0)).array().clipUnder(0.0);
    }
  } // for each row
}

void Meddis::reset(){
  // Initial values
  LOUDIA_DEBUG("MEDDIS: Resetting...");

  kt = MatrixXR::Constant(1, _channels, g * A / (A + B));
  spont = kt * (M * y);
  spont = spont.cwise() / ( (kt*l).cwise() + (y * ( l + r )) );

  c = spont;
  q = (c * ( l + r )).cwise() / kt;
  w = c * r / x;
}


int Meddis::channels() const {
  return _channels;
}

Real Meddis::sampleRate() const {
  return _sampleRate;
}
