from __future__ import absolute_import

from pbcore.io import BasH5Reader, BasH5Collection

__author__ = 'Brett Bowman'

def loadFromFile(filename):
    try:
        reader = BasH5Collection(filename)
    except:
        raise IOError("Invalid Bas.H5 file or collection")
    return reader

def getWhiteListType(filename):
    if filename.endswith('.csv'):
        return "CSV"
    elif filename.endswith('.txt', '.zmw'):
        return "ZMW"
    else:
        raise TypeError("Invalid Whitelist filetype")

def readZmwWhiteList(filename):
    try:
        zmws = set([l.strip() for l in open(filename) if len(l.strip()) > 1])
    except:
        raise IOError("Invalid Whitelist file")
    return zmws

def readCsvWhiteList(filename):
    try:
        zmws = set([l.strip().split(',')[0] for l in open(filename) if len(l.strip()) > 1])
    except:
        raise IOError("Invalid Whitelist file")
    return zmws

def readWhiteListIdx(filename):
    try:
        idx = {l.strip().split(',')[0]:l.strip().split(',')[1] for l in open(filename) if len(l.strip()) > 1}
    except:
        raise IOError("Invalid Whitelist file")
    return idx

def filterZmws(collection, whiteList):
    return [collection[z] for z in whiteList]
