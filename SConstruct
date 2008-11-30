"""                                                        
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
"""


VariantDir('build', 'src')

env = Environment()  # Initialize the environment

# Add the includes
env.Append(CPPPATH = ['/home/rmarxer/dev/eigen2', './build'])

# Add the lib dir
env.Append(LIBPATH = ['./build'])

# Add the flags for the compiler
env.Append(CXXFLAGS = '-O3 -DRICAUDIO_DEBUG')
#env.Append(CXXFLAGS = '-O3 -msse2 -DEIGEN_NO_DEBUG -DRICAUDIO_DEBUG')

# Build the tests
env.Program('build/test_meddis', ['build/tests/test_meddis.cpp', 'build/meddis.cpp', 'build/debug.cpp'])
env.Program('build/test_filter', ['build/tests/test_filter.cpp', 'build/filter.cpp', 'build/debug.cpp'])
env.Program('build/test_spectralbands', ['build/tests/test_spectralbands.cpp', 'build/spectralbands.cpp', 'build/debug.cpp'])
env.Program('build/test_melbands', ['build/tests/test_melbands.cpp', 'build/melbands.cpp', 'build/spectralbands.cpp', 'build/debug.cpp'])
env.Program('build/test_dct', ['build/tests/test_dct.cpp', 'build/dct.cpp', 'build/debug.cpp'])


env['SHLIBPREFIX'] = ''
env.Append(CPPPATH = [ '/usr/include/python2.5',
                       '/usr/lib/python2.5/site-packages/numpy/core/include/numpy/' ])
env.SharedLibrary('build/pycricaudio', ['build/python/cricaudio.cpp',
                                        'build/python/pymeddis.cpp', 'build/meddis.cpp',
                                        'build/python/pyfilter.cpp', 'build/filter.cpp',
                                        'build/debug.cpp'])

