#! /usr/bin/env python
# encoding: utf-8
# Ricard Marxer 2008

VERSION='0.0.1'
APPNAME='ricaudio'

# these variables are mandatory,
srcdir = '.'
blddir = 'build'

def set_options(opt):
	pass

def configure(conf):
	pass

def build(bld):
        bld.add_subdirs('src')
        
        
