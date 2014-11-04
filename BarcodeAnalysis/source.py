from __future__ import absolute_import

from pbcore.io import BasH5Reader, BasH5Collection

__author__ = 'Brett Bowman'

def loadFromFile(filename):
    try:
        reader = BasH5Collection(filename)
    except:
        raise IOError("Invalid Bas.H5 file or collection")
    return reader

def readWhiteList(filename):
    try:
        zmws = set([l.strip() for l in open(filename) if len(l.strip()) > 1])
    except:
        raise IOError("Invalid Whitelist file")
    return zmws

def filterZmws(collection, whiteList):
    return [z for z in collection if z.zmwName in whiteList]
