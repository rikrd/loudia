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
    def __init__(self, filename, firstChar = 10, lastChar = 100):
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
        self.firstChar = firstChar
        self.lastChar = lastChar
        
    def analyze(self):
        for code in range(self.firstChar, self.lastChar):
            char = Character(code)
            
            glyph = Glyph(char, self.font)
            glyph.analyze()
            
            self.chars[char] = glyph.features

class Glyph:
    def __init__(self, char, font, numPoints = 2):
        self.font = font
        self.char = char
        self.features = {}
        self.numPoints = numPoints

    def analyze(self):
        self.analyze_metrics()
        self.analyze_shape()
        
    def analyze_metrics(self):
        # Create the metrics object
        self.metrics = QFontMetricsF( self.font )

        # Extract features
        self.features['height'] = self.metrics.height()
        self.features['width'] = self.metrics.width(self.char)
        self.features['leading'] = self.metrics.leading()
        return

    def analyze_shape(self):
        # Create the shape object
        self.shape = QPainterPath()
        self.shape.addText(0.0, 0.0, self.font, QString(self.char))
        
        # Extract features
        self.features['angles'] = [self.shape.angleAtPercent(t/float(self.numPoints)) for t in range(self.numPoints)]
        return

        
class Character(QChar):
    def __init__(self, code):
        QChar.__init__( self, code )

    def direction(self):
        return CharDirection()[QChar.direction(self)]

    def category(self):
        return CharCategory()[QChar.category(self)]        
        
if __name__ == '__main__':    
    app = QApplication( sys.argv )

    charCategories = CharCategory()
    charDirections = CharDirection()
    print charCategories
    print charDirections
    
    fe = Font( sys.argv[1] )
    fe.analyze()
    print [(str(char.toLatin1()), char.direction(), char.category(), value) for char, value in fe.chars.iteritems()]
    #app.exec_()
