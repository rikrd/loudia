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

import distutils.sysconfig

VariantDir('build', 'src')

env = Environment()  # Initialize the environment

# Add the includes
env.Append(CPPPATH = ['/home/rmarxer/dev/eigen2', './build'])

# Add the lib dir
env.Append(LIBPATH = ['./build'])

# Add the flags for the compiler
#env.Append(CXXFLAGS = '-O3 -DRICAUDIO_DEBUG')
env.Append(CXXFLAGS = '-O3 -fPIC -msse2 -DRICAUDIO_DEBUG')
#env.Append(CXXFLAGS = '-O3 -msse2 -DEIGEN_NO_DEBUG')

env.Append(LINKFLAGS = '-lfftw3f')

# Build the tests
env.Program('build/test_meddis', ['build/tests/test_meddis.cpp', 'build/meddis.cpp', 'build/debug.cpp'])
env.Program('build/test_filter', ['build/tests/test_filter.cpp', 'build/filter.cpp', 'build/debug.cpp'])
env.Program('build/test_spectralbands', ['build/tests/test_spectralbands.cpp', 'build/spectralbands.cpp', 'build/debug.cpp'])
env.Program('build/test_melbands', ['build/tests/test_melbands.cpp', 'build/melbands.cpp', 'build/spectralbands.cpp', 'build/debug.cpp'])
env.Program('build/test_dct', ['build/tests/test_dct.cpp', 'build/dct.cpp', 'build/debug.cpp'])
env.Program('build/test_mfcc', ['build/tests/test_mfcc.cpp', 'build/mfcc.cpp', 'build/debug.cpp', 'build/melbands.cpp', 'build/spectralbands.cpp', 'build/dct.cpp'])
env.Program('build/test_aok', ['build/tests/test_aok.cpp', 'build/aok.cpp', 'build/debug.cpp'])
env.Program('build/test_fft', ['build/tests/test_fft.cpp', 'build/fft.cpp', 'build/debug.cpp'])
env.Program('build/test_window', ['build/tests/test_window.cpp', 'build/window.cpp', 'build/debug.cpp'])
env.Program('build/test_spectralreassignment', ['build/tests/test_spectralreassignment.cpp', 'build/spectralreassignment.cpp', 'build/window.cpp', 'build/fft.cpp', 'build/debug.cpp'])


# Build python bindings
env.Append(CPPPATH = [distutils.sysconfig.get_python_inc(), '/usr/lib/python2.5/site-packages/numpy/core/include/numpy/'],
           SHLIBPREFIX="")

# Build manual bindings
env.SharedLibrary('build/ricaudio', ['build/python/ricaudio.cpp',
                                     'build/python/pymeddis.cpp', 'build/meddis.cpp',
                                     'build/python/pymfcc.cpp', 'build/mfcc.cpp', 'build/melbands.cpp',
                                     'build/spectralbands.cpp', 'build/dct.cpp',
                                     'build/python/pyaok.cpp', 'build/aok.cpp',
                                     'build/python/pyfilter.cpp', 'build/filter.cpp',
                                     'build/python/pyfft.cpp', 'build/fft.cpp',
                                     'build/python/pywindow.cpp', 'build/window.cpp',
                                     'build/python/pyspectralreassignment.cpp', 'build/spectralreassignment.cpp',
                                     'build/debug.cpp'], SHLIBPREFIX="")


# Build SWIG generated python bindings
"""
env.Append(SWIGFLAGS=['-c++', '-python', '-Wall', '-I./build', '-I/home/rmarxer/dev/eigen2'])
env.SharedLibrary('build/swig/_ricaudio.so', ['build/swig/ricaudio_wrap.cc', 'build/swig/ricaudio.i'], SHLIBPREFIX="")
"""
