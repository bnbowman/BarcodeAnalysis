from __future__ import absolute_import

__author__ = 'Brett Bowman'

import string
import h5py
import numpy as np

from pbcore.io import FastaRecord
import ConsensusCore as cc

fst = lambda t: t[0]
snd = lambda t: t[1]

QUIVER_FEATURES = ["InsertionQV",
                   "SubstitutionQV",
                   "DeletionQV",
                   "DeletionTag",
                   "MergeQV"]

REV_COM = string.maketrans("AGCT-", "TCGA-")

def reverse_complement( sequence ):
    if isinstance( sequence, FastaRecord ):
        return FastaRecord( sequence.name,
                            rev_com_str(sequence.sequence) )
    if isinstance( sequence, str ):
        return rev_com_str( sequence )

def rev_com_str( seq ):
    return seq.translate( REV_COM )[::-1]

def asFloatFeature( array ):
    return cc.FloatFeature(np.array(array, dtype=np.float32))

def arrayFromDataset(ds, offsetBegin, offsetEnd):
    """
    Extract a one-dimensional array from an HDF5 dataset.
    """
    shape = (offsetEnd - offsetBegin,)
    a = np.ndarray(shape=shape, dtype=ds.dtype)
    mspace = h5py.h5s.create_simple(shape)
    fspace = ds.id.get_space()
    fspace.select_hyperslab((offsetBegin,), shape, (1,))
    ds.id.read(mspace, fspace, a)
    return a

def die(msg):
    print msg
    sys.exit(-1)
