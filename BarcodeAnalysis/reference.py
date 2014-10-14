from __future__ import absolute_import

from pbcore.io import FastaReader

__author__ = 'Brett Bowman'

def loadFromFile(filename):
    """
    Reads barcode sequences from a FASTA file into a dictionary so
    that individual records could be accessed quickly
    """
    try:
        reader = FastaReader(filename)
    except:
        raise IOError("Invalid FASTA file")

    return {r.name.strip().split()[0]: r for r in reader}
