/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#include "meddis.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <iostream>

using namespace std;

int main() {
  Meddis mds(8000, 30);
  mds.setup();

  int rows = 1024;
  int cols = 30;
  MatrixXR result(rows, cols);

  MatrixXR in(rows, cols);
  
  for (int j = 0; j < rows; j++) {
    for (int k = 0; k < cols; k++) {
      in(j, k) = 3.0;
    }
  }
  
  for (int i=0; i<1000; i++) {   
    mds.process(in, &result);
  }

  return 0;
}

