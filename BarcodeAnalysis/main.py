#!/usr/bin/env python
from __future__ import absolute_import

__author__ = "Brett Bowman"

import logging, string, math
from collections import defaultdict
from pkg_resources import Requirement, resource_filename
from pbcore.io import FastaReader

import h5py
import numpy as np

from BarcodeAnalysis._version import __version__
from BarcodeAnalysis import reference, source, configure
from BarcodeAnalysis.options import (options,
                                     parseOptions,
                                     consensusCoreVersion)
from BarcodeAnalysis.ConsensusCoreRead import ConsensusCoreRead
from BarcodeAnalysis.utils import QUIVER_FEATURES, reverse_complement

import ConsensusCore as cc

def _getAbsPath(filename):
    return resource_filename(Requirement.parse('BarcodeAnalysis'), 'BarcodeAnalysis/{0}'.format( filename ))

def _getDefaultParametersFile():
    return _getAbsPath('QuiverParameters.ini')

def barcodeCsvLine( holeNum, scores, adapters ):
    sortedScores = sorted(scores.iteritems(), key=lambda x: x[1], reverse=True)
    bestIdx, bestScore = sortedScores[0]
    secondIdx, secondScore = sortedScores[1] if len(sortedScores) > 1 else ("N/A", 0)
    bestProb = math.exp( bestScore )
    secondProb = math.exp( secondScore ) if secondScore else 1e-100
    ratio = bestProb/secondProb
    return "{0},{1},{2},{3},{4},{5},{6},{7}\n".format(holeNum, adapters,
                                                      bestIdx, bestScore, bestProb,
                                                      secondIdx, secondScore, secondProb,
                                                      ratio)

class BarcodeAnalyzer(object):
    """
    A tool for analyzing barcode sequences with Quiver
    """
    def __init__(self):
        self._inputReader = None
        self._sequenceZmws = None
        self._barcodeDict = None
        self._barcodeSequences = None
        self._barcodeNames = None
        self._barcodePairs = None
        self._configDict = None
        self._windowSize = 25
        self._adapterPad = 0

    def _setupLogging(self):
        if options.quiet:
            logLevel = logging.ERROR
        elif options.verbosity >= 2:
            logLevel = logging.DEBUG
        elif options.verbosity == 1:
            logLevel = logging.INFO
        else:
            logLevel = logging.WARNING
        logFormat = '[%(levelname)s] %(message)s'
        logging.basicConfig(level=logLevel, format=logFormat)

    def _loadData(self):
        logging.info("Loading input data")
        reader = source.loadFromFile(options.inputFilename)
        if options.nZmws < 0:
            self._sequencingZmws = reader.sequencingZmws
        else:
            self._sequencingZmws = reader.sequencingZmws[:options.nZmws]
        self.inputReader = reader

    def _loadBarcodes(self):
        """
        Read barcode names, sequences, and pairs from the barcode Fasta file
        """
        logging.info("Loading barcodes")
        raw_seqs = set()    # 1. Make sure no barcode sequences are duplicate
        names = []          # 2. Keep a list of all unique barcode names
        sequences = {}      # 3. Create a dictionary of all barcode sequences
        for barcode in FastaReader(options.barcodeFilename):
            name = barcode.name.strip().split()[0]
            # Check the barcode name
            if name in names:
                raise ValueError("Duplicate barcode name in '%s'".format(name))
            else:
                names.append( name )

            # Check the forward seqeunce
            if barcode.sequence in raw_seqs:
                raise ValueError("Duplicate barcode sequence in '%s'".format(name))
            else:
                raw_seqs.add( barcode.sequence )

            # Check the reverse complement sequence
            rc_barcode = reverse_complement( barcode )
            if rc_barcode.sequence in raw_seqs:
                raise ValueError("Duplicate barcode sequence in '%s'".format(name))
            else:
                raw_seqs.add( rc_barcode.sequence )

            # If both pass, add the sequences and pair-wise combinations
            sequences[(name, 'FORWARD')] = barcode.sequence
            sequences[(name, 'REVERSE')] = rc_barcode.sequence
        self._barcodeSequences = sequences
        self._barcodeNames = set(names)
        self._barcodePairs = [(names[i], names[i+1]) for i in range(0,len(names)-1,2)]

    @property
    def windowSize(self):
        return self._windowSize
    @property
    def adapterPad(self):
        return self._adapterPad
    @property
    def chemistry(self):
        return self.inputReader.sequencingChemistry
    @property
    def model(self):
        return "AllQVsMergingByChannelModel"
    @property
    def quiverModel(self):
        return "{0}.{1}".format(self.chemistry, self.model)
    @property
    def parametersFile(self):
        if options.parametersFile is None:
            return _getDefaultParametersFile()
        else:
            return options.parametersFile
    @property
    def quiverConfig(self):
        if self._configDict is None:
            self._configDict = configure.loadQuiverConfigs(self.parametersFile)
        try:
            config = self._configDict[self.quiverModel]
        except:
            raise KeyError("No QuiverConfig found for '{0}'".format(self.quiverModel))
        return config
    @property
    def readScorer(self):
        return cc.ReadScorer(self.quiverConfig)
    @property
    def barcodeSequences(self):
        # If the sequences haven't been set, calculate them
        if self._barcodeSequences is None:
            self._prepareBarcodeSequences()
        # Return an iterator over those sequences
        return self._barcodeSequences

    def _getBaxForHole( self, holeNum ):
        partNum = self.inputReader._holeLookup(holeNum) - 1
        return self.inputReader._parts[partNum]

    def _makeRead( self, holeNum, start, end ):
        bax = self._getBaxForHole( holeNum )
        return ConsensusCoreRead(bax, holeNum, start, end, self.chemistry)

    def scoreBarcode(self, barcode, window):
        return self.readScorer.Score(barcode, window.read)

    def scoreWindow(self, window):
        scores = {}
        for key, barcode in self.barcodeSequences.iteritems():
            scores[key] = self.scoreBarcode(barcode, window)
        return scores

    def scoreLeftWindow( self, zmw, adapter ):
        hqStart, hqEnd = zmw.hqRegion
        adpStart, adpEnd = adapter
        windowStart = max(hqStart, adpStart - self.windowSize)
        windowEnd = min(hqEnd, adpStart + self.adapterPad)
        window = self._makeRead( zmw.holeNumber, windowStart, windowEnd )
        return self.scoreWindow( window )

    def scoreRightWindow( self, zmw, adapter ):
        hqStart, hqEnd = zmw.hqRegion
        adpStart, adpEnd = adapter
        windowStart = max(hqStart, adpEnd - self.adapterPad)
        windowEnd = min(hqEnd, adpEnd + self.windowSize)
        window = self._makeRead( zmw.holeNumber, windowStart, windowEnd )
        return self.scoreWindow( window )

    def combineScores( self, left, right ):
        if left is None:
            return right
        elif right is None:
            return left
        else:
            combined = {}
            for name in self._barcodeNames:
                forwardTup = (name, 'FORWARD')
                reverseTup = (name, 'REVERSE')
                combined[forwardTup] = left[forwardTup] + right[reverseTup]
                combined[reverseTup] = right[forwardTup] + left[reverseTup]
            return combined

    def bestScoreByBarcode( self, scores ):
        bestScores = {}
        for name in self._barcodeNames:
            forwardScore = scores[(name, 'FORWARD')]
            reverseScore = scores[(name, 'REVERSE')]
            if forwardScore > reverseScore:
                bestScores[name] = forwardScore
            else:
                bestScores[name] = reverseScore
        return bestScores

    def addScores( self, oldScores, newScores ):
        for key, score in newScores.iteritems():
            oldScores[key].append( score )
        return oldScores

    def scorePairs( self, scores ):
        pairScores = {}
        for forward, reverse in self._barcodePairs:
            pairName = '%s--%s' % (forward, reverse)
            pairScores[pairName] = np.maximum(scores[forward], scores[reverse])
        return pairScores

    def sumScores( self, scores ):
        finalScores = {k: sum(v) for k, v in scores.iteritems()}
        return finalScores

    def scoreZmw(self, zmw):
        """Combine and report the score from each adapter in the ZMW
        """
        # Iterate over each adapter, scoring them separately
        totalScores = defaultdict(list)
        for adapter in zmw.adapterRegions:
            # Score the left and right windows independently and combine them
            leftScores  = self.scoreLeftWindow( zmw, adapter )
            rightScores = self.scoreRightWindow( zmw, adapter )
            combinedScores = self.combineScores( leftScores, rightScores )
            # Take the best score for each barcode from the two combined options
            bestScores = self.bestScoreByBarcode( combinedScores )
            totalScores = self.addScores( totalScores, bestScores)
        # Take the best score of each barcode pair for each position
        pairScores = self.scorePairs( totalScores )
        finalScores = self.sumScores( pairScores )
        # Return the combined scores
        return finalScores

    def writeZmwScores(self, zmwScores, zmwAdapters):
        with open(options.outputFilename, 'w') as handle:
            handle.write("HoleNumber,NumAdapters,Idx1,LogScore1,Psudo1,Idx2,LogScore2,Pseudo2,Ratio\n")
            for holeNum in sorted(zmwScores.keys()):
                barcodeScores = zmwScores[holeNum]
                adapterCount = zmwAdapters[holeNum]
                handle.write( barcodeCsvLine(holeNum, barcodeScores, adapterCount) )

    def main(self):
        parseOptions()
        self._setupLogging()

        logging.info("h5py version: %s" % h5py.version.version)
        logging.info("hdf5 version: %s" % h5py.version.hdf5_version)
        logging.info("ConsensusCore version: %s" % consensusCoreVersion())
        logging.info("BarcodeAnalysis version: %s" % __version__)

        logging.info("Starting.")
        self._loadData()
        self._loadBarcodes()

        zmwScores = {}
        zmwAdapters = {}
        for holeNum in self._sequencingZmws:
            zmw = self.inputReader[holeNum]
            if len(zmw.adapterRegions) == 0:
                logging.debug("Skipping ZMW #{0} -- no adapters".format(holeNum))
                continue
            logging.debug("Labelling ZMW #{0}".format(holeNum))
            scores = self.scoreZmw( zmw )
            zmwScores[zmw.holeNumber] = scores
            zmwAdapters[zmw.holeNumber] = len(zmw.adapterRegions)

        self.writeZmwScores( zmwScores, zmwAdapters )

        return 0

    def printScores(self, scores):
        for key, value in sorted(scores.iteritems()):
            print key, value, math.exp( value )
        print

    def printScoreArrays(self, scores):
        for key, values in scores.iteritems():
            print key, values, math.exp( sum(values) )
        print

def main():
    return BarcodeAnalyzer().main()
