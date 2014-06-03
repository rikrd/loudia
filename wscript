#! /usr/bin/env python
# encoding: utf-8
# Ricard Marxer 2008

# the following two variables are used by the target "waf dist"
APPNAME = 'loudia'
VERSION = '0.1-dev'

# these variables are mandatory,
srcdir = '.'
blddir = 'build'

def options(opt):
        opt.recurse('src')
        opt.recurse('doc')

        opt.add_option('--debug', action='store_true', default=False, help='Compile in debug mode')
        opt.add_option('--no-python-bindings', action='store_true', default=False, help='Generate and compile the Python bindings')
        opt.add_option('--cpptests', action='store_true', default=False, help='Compile C++ tests.')
        opt.add_option('--doc', action='store_true', default=False, help='Generate the documentation')

def configure(conf):
        from waflib import Options
        
        conf.env['option_doc'] = Options.options.doc
        conf.env['option_debug'] = Options.options.debug
        conf.env['option_no_python_bindings'] = Options.options.no_python_bindings
        conf.env['option_cpptests'] = Options.options.cpptests

        conf.recurse('src')

        if conf.env['option_doc']:
                conf.recurse('doc')


def build(bld):
        bld.recurse('src')

        if bld.env['option_doc']:
                bld.recurse('doc')

