#! /usr/bin/env python
# encoding: utf-8
# Ricard Marxer 2008

# the following two variables are used by the target "waf dist"
APPNAME = 'ricaudio'
VERSION = '0.1-dev'

# these variables are mandatory,
srcdir = '.'
blddir = 'build'

def set_options(opt):
        opt.sub_options('src')
        #opt.sub_options('doc')

def configure(conf):
        conf.sub_config('src')
        #conf.sub_config('doc')
        
def build(bld):
        bld.add_subdirs('src')
        #bld.add_subdirs('doc')        
        
