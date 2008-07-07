#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import sys

class CharCategory(dict):
    def __init__(self):
        # Create a dictionary containing:
        #   key:    the value of the QChar.Category
        #   value:  category hierarchy in a list form
        self.update(dict([(getattr(QChar, key), key.split('_')) for key in dir(QChar) if type(getattr(QChar, key)) == type(QChar.Symbol_Currency)]))

class CharDirection(dict):
    def __init__(self):
        # Create a dictionary containing:
        #   key:    the value of the QChar.Category
        #   value:  category hierarchy in a list form
        self.update(dict([(getattr(QChar, key), key.split('_')) for key in dir(QChar) if type(getattr(QChar, key)) == type(QChar.DirAL)]))

class Font:
    def __init__(self, filename):
        self.features = {}
        
        self.filename = filename
        self.db = QFontDatabase()
        before = set( self.db.families() )
        self.db.removeAllApplicationFonts()
        self.db.addApplicationFont( filename )
        after = set( self.db.families() )
        self.family = list( after - before )[0]

        print 'Loading font "%s" ...' % str( self.family )
        
        self.font = QFont( self.family )
        self.chars = {}
        self.numChars = 10
        
    def analyze(self):
        for code in range(self.numChars):
            char = Character(code)
            
            glyph = Glyph(char, self.font)
            glyph.analyze()
            
            self.chars[char] = glyph.features

class Glyph:
    def __init__(self, char, font):
        self.font = font
        self.char = char
        self.features = {}

    def analyze(self):
        self.analyze_metrics()
        self.analyze_shape()
        
    def analyze_metrics(self):
        self.metrics = QFontMetricsF( self.font )
        
        self.features['height'] = self.metrics.height()
        self.features['width'] = self.metrics.width(self.char)
        self.features['leading'] = self.metrics.leading()
        return

    def analyze_shape(self):
        self.shape = QPainterPath()
        self.shape.addText(0.0, 0.0, self.font, QString(self.char))
        print self.shape.angleAtPercent(0.5)
        return

        
class Character(QChar):
    def __init__(self, code):
        QChar.__init__( self, code )
        
if __name__ == '__main__':    
    app = QApplication( sys.argv )

    charCategories = CharCategory()
    charDirections = CharDirection()
    print charCategories
    print charDirections

    '''
    for i in range( 1000 ):
        ch = Character( i )
        print charCategories[ch.category()]
        print charDirections[ch.direction()]
    '''
    
    fe = Font( sys.argv[1] )
    fe.analyze()
    print [(str(key.char.toLatin1()), value) for key, value in fe.chars.iteritems()]
    #app.exec_()
