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

#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.h"
#define EIGEN_CWISE_PLUGIN "CwiseAddons.h"
#define EIGEN_FUNCTORS_PLUGIN "FunctorsAddons.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

using namespace std;


int main() {
  Eigen::MatrixXf m = Eigen::MatrixXf::Random(4, 4);
  Eigen::MatrixXf n = Eigen::MatrixXf::Random(4, 4);
 
  cout << "-----------Test 1-------------" << endl;

  cout << "A random matrix: " << endl;

  cout << m << endl;
  
  cout << "The maxed matrix:" << endl;

  cout << m.cwise().clipUnder(0.4) << endl;

  cout << "-----------Test 2-------------" << endl;

  cout << "The first row: " << endl;

  cout << m.row(0) << endl;
  
  cout << "The maxed first row:" << endl;

  cout << m.row(0).cwise().clipUnder(0.4) << endl;

  cout << "-----------Test 3-------------" << endl;

  cout << "The first row, binary op: " << endl;

  cout << m.row(0) + n.row(0) << endl;
  
  cout << "The maxed first row, binary op:" << endl;

  cout << (m.row(0) + n.row(0)).cwise().clipUnder(0.4) << endl;

}
