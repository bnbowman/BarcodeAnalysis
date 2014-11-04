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

def barcodeCsvLine( zmw, scores ):
    adapters = len(zmw.adapterRegions)
    sortedScores = sorted(scores.iteritems(), key=lambda x: x[1], reverse=True)
    bestIdx, bestScore = sortedScores[0]
    secondIdx, secondScore = sortedScores[1] if len(sortedScores) > 1 else ("N/A", 0)
    bestProb = math.exp( bestScore )
    secondProb = math.exp( secondScore ) if secondScore else 1e-100
    try:
        ratio = bestProb/secondProb
    except ZeroDivisionError:
        ratio = 1e-100
    return "{0},{1},{2},{3},{4},{5},{6},{7}\n".format(zmw.zmwName, adapters,
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
        self._barcodeLength = None
        self._whitelist = None
        self._configDict = None
        self._chemistry = None
        self._windowSize = 25
        self._adapterPad = 0
        self._startTimeCutoff = 10.0

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

    def _loadChemistry(self):
        all_chemistries = [self.inputReader[k].sequencingChemistry for k in self.inputReader.readers.keys()]
        used_chemistries = sorted(set(all_chemistries))
        print all_chemistries
        print used_chemistries
        assert len(used_chemistries) == 1
        self._chemistry = used_chemistries[0]

    def _loadData(self):
        logging.info("Loading input data")
        reader = source.loadFromFile(options.inputFilename)
        # Filter out only the whitelisted ZMWs if a whitelist was given
        print "opened"
        if options.whiteList:
            whitelist = source.readWhiteList(options.whiteList)
            print "whitelist"
            sequencingZmws = source.filterZmws(reader, whitelist)
            print "subset"
        else:
            sequencingZmws = [z.zmwName for z in reader]

        # Clip the list of ZMWs if nZMWs was set
        if options.nZmws < 0:
            self._sequencingZmws = sequencingZmws
        else:
            self._sequencingZmws = sequencingZmws[:options.nZmws]

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
        # Verify that all of the barcodes are the same length
        bc_lengths = list(set([len(s) for s in sequences.itervalues()]))
        if len(bc_lengths) > 1:
            msg = "Multiple barcode lengths detected - {0}".format(bc_lengths)
            logging.error( msg )
            raise ValueError( msg )
        self._barcodeLength = bc_lengths[0]
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
        return self._chemistry
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

    def _makeRead( self, zmw, start, end ):
        if (end-start) < self._barcodeLength:
            logging.debug("ZMW #{0} - Skipping barcode window ({1}, {2}) - too small".format(zmw.holeNumber, start, end))
            return None
        #bax = self._getBaxForHole( holeNum )
        return ConsensusCoreRead(zmw.baxH5, zmw.holeNumber, start, end, self.chemistry)

    def scoreBarcode(self, barcode, window):
        return self.readScorer.Score(barcode, window.read)

    def scoreWindow(self, zmw, start, end):
        window = self._makeRead( zmw, start, end )
        # If there isn't enough data to make a read from, return None
        if window is None:
            return None
        # Otherwise score the window against all possible barcodes
        else:
            scores = {}
            for key, barcode in self.barcodeSequences.iteritems():
                scores[key] = self.scoreBarcode(barcode, window)
            return scores

    def scoreLeftWindow( self, zmw, adapter ):
        hqStart, hqEnd = zmw.hqRegion
        adpStart, adpEnd = adapter
        windowStart = max(hqStart, adpStart - self.windowSize)
        windowEnd = min(hqEnd, adpStart + self.adapterPad)
        #windowStart = max(hqStart, adpStart + 20)
        #windowEnd = min(hqEnd, adpStart + 45)
        return self.scoreWindow( zmw, windowStart, windowEnd )

    def scoreRightWindow( self, zmw, adapter ):
        hqStart, hqEnd = zmw.hqRegion
        adpStart, adpEnd = adapter
        windowStart = max(hqStart, adpEnd - self.adapterPad)
        windowEnd = min(hqEnd, adpEnd + self.windowSize)
        #windowStart = max(hqStart, adpEnd - 20)
        #windowEnd = min(hqEnd, adpEnd - 45)
        return self.scoreWindow( zmw, windowStart, windowEnd )

    def _getWindowReads( self, zmw ):
        hqStart, hqEnd = zmw.hqRegion
        for adapter in zmw.adapterRegions:
            adpStart, adpEnd = adapter
            leftStart = adpStart - self.windowSize
            leftEnd = adpStart + self.adapterPad
            if leftStart >= hqStart:
                yield self._makeRead( zmw, leftStart, leftEnd )
            else:
                yield None
            rightStart = adpEnd - self.adapterPad
            rightEnd = adpEnd + self.windowSize
            if rightEnd <= hqEnd:
                yield self._makeRead( zmw, rightStart, rightEnd )
            else:
                yield None

    def scoreFirstWindow(self, zmw):
        s = zmw.zmwMetric('HQRegionStartTime')
        e = zmw.zmwMetric('HQRegionEndTime')
        # s<e => has HQ.
        if s < e and s <= self._startTimeCutoff:
            end = self.windowSize + self.adapterPad
            end = end if zmw.hqRegion[1] > end else zmw.hqRegion[1]
            scores = self.scoreWindow( zmw.holeNumber, 0, end )
            bestScores = self.bestScoreByBarcode( scores )
            pairScores = self.scorePairs( bestScores )
        else:
            logging.debug("Skipping ZMW #{0} -- start-time too late".format(zmw.holeNumber))
            pairScores = None
        return pairScores

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
            if combinedScores is None:
                continue
            # Take the best score for each barcode from the two combined options
            bestScores = self.bestScoreByBarcode( combinedScores )
            totalScores = self.addScores( totalScores, bestScores)
        # Take the best score of each barcode pair for each position
        pairScores = self.scorePairs( totalScores )
        finalScores = self.sumScores( pairScores )
        # Return the combined scores
        return finalScores

    def openOutputFile(self):
        if not options.appendOutput:
            handle = open(options.outputFilename, 'w')
            handle.write("Zmw,NumAdapters,Idx1,LogScore1,Psudo1,Idx2,LogScore2,Pseudo2,Ratio\n")
        else:
            handle = open(options.outputFilename, 'a')
        return handle

    def main(self):
        parseOptions()
        self._setupLogging()

        logging.info("h5py version: %s" % h5py.version.version)
        logging.info("hdf5 version: %s" % h5py.version.hdf5_version)
        logging.info("ConsensusCore version: %s" % consensusCoreVersion())
        logging.info("BarcodeAnalysis version: %s" % __version__)

        logging.info("Starting.")
        self._loadData()
        self._loadChemistry()
        self._loadBarcodes()

        if options.tSNE:
            with self.openOutputFile() as handle:
                for holeNum in self._sequencingZmws:
                    zmw = self.inputReader[holeNum]
                    if len(zmw.adapterRegions) >= 3:
                        for window in self._getWindowReads(zmw):
                            handle.write( window.to_csv + '\n' )
        elif options.adapterSizes:
            with self.openOutputFile() as handle:
                for zmw in self._sequencingZmws:
                    lengths = [(len(w) if w else 0) for w in self._getWindowReads(zmw)]
                    print '%s,%s,%s' % (zmw.zmwName, ','.join([str(l) for l in lengths]), min(lengths))
        else:
            with self.openOutputFile() as handle:
                for zmw in self._sequencingZmws:
                    # If scoreFirst is on, ONLY score the first
                    if options.scoreFirst and len(zmw.adapterRegions) == 0:
                        logging.debug("Labelling 5'-end of ZMW #{0} -- no adapters".format(zmw.holeNumber))
                        scores = self.scoreFirstWindow( zmw )
                    # Otherwise score normally
                    else:
                        if len(zmw.adapterRegions) == 0:
                            logging.debug("Skipping ZMW #{0} -- no adapters".format(zmw.holeNumber))
                            continue
                        logging.debug("Labelling ZMW #{0}".format(zmw.holeNumber))
                        scores = self.scoreZmw( zmw )
                    print scores
                    if scores is not None:
                        handle.write( barcodeCsvLine( zmw, scores ) )

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
