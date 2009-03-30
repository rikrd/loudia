#! /usr/bin/env python
# encoding: utf-8
# Ricard Marxer 2008

# the following two variables are used by the target "waf dist"
APPNAME = 'loudia'
VERSION = '0.1-dev'

# these variables are mandatory,
srcdir = '.'
blddir = 'build'

def set_options(opt):
        opt.sub_options('src')
        opt.sub_options('doc')

        opt.add_option('--debug', action='store_true', default=False, help='Compile in debug mode')
        opt.add_option('--no-python-bindings', action='store_true', default=False, help='Generate and compile the Python bindings')
        opt.add_option('--cpptests', action='store_true', default=False, help='Compile C++ tests.')
	opt.add_option('--doc', action='store_true', default=False, help='Generate the documentation')

def configure(conf):
        conf.sub_config('src')

        import Options
        if Options.options.doc:
                conf.sub_config('doc')
        
def build(bld):
        bld.add_subdirs('src')

        import Options
        if Options.options.doc:
                bld.add_subdirs('doc')
        
