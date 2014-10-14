from __future__ import absolute_import

from pbcore.io import BasH5Reader, BasH5Collection

__author__ = 'Brett Bowman'

def loadFromFile(filename):
    """
    Reads barcode sequences from a FASTA file into a dictionary so
    that individual records could be accessed quickly
    """
    try:
        reader = BasH5Reader(filename)
    except:
        raise IOError("Invalid Bas.H5 file")

    return reader
