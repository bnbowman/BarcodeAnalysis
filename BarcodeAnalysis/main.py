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
#from BarcodeAnalysis.ConsensusCoreRead import ConsensusCoreRead
from BarcodeAnalysis.utils import QUIVER_FEATURES, reverse_complement
from BarcodeAnalysis.BarcodeScorer import BarcodeScorer

#import ConsensusCore as cc

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

def barcodeCsvLine2( zmw, trueIdx, scores ):
    adapters = len(zmw.adapterRegions)
    sortedScores = sorted(scores.iteritems(), key=lambda x: x[1], reverse=True)
    bestIdx, bestScore = sortedScores[0]
    secondIdx, secondScore = sortedScores[1] if len(sortedScores) > 1 else ("N/A", 0)
    bestAvg = bestScore / float(adapters)
    secondAvg = secondScore / float(adapters)
    ratio = bestScore/float(secondScore)
    bestTrue = 1 if bestIdx == trueIdx else 0
    secondTrue = 1 if secondIdx == trueIdx else 0
    return "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".format(zmw.zmwName, trueIdx, adapters,
                                                                        bestIdx, bestScore, bestAvg, bestTrue,
                                                                        secondIdx, secondScore, secondAvg, secondTrue,
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
        self._barcodePairNames = None
        self._barcodeLength = None
        self._whiteList = None
        self._whiteListIdx = None
        self._configDict = None
        self._chemistry = None
        self._windowSize = 25
        self._adapterPad = 0
        self._insertPad = 9
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

    def _loadWhiteList(self):
        logging.info("Loading whitelist information")
        if options.whiteList:
            whiteListType = source.getWhiteListType( options.whiteList )
            if whiteListType == "ZMW":
                self._whiteList = source.readZmwWhiteList( options.whiteList )
            elif whiteListType == "CSV":
                self._whiteList = source.readCsvWhiteList( options.whiteList )
                self._whiteListIdx = source.readWhiteListIdx( options.whiteList )
            elif whiteListType == "FASTA":
                self._whiteList = source.readFastaWhiteList( options.whiteList )

    def _loadData(self):
        logging.info("Loading input data")
        self._inputReader = source.loadFromFile(options.inputFilename)

        # Filter out only the whitelisted ZMWs if a whitelist was given
        if options.whiteList:
            sequencingZmws = source.filterZmws(self._inputReader, self._whiteList)
        else:
            sequencingZmws = [z.zmwName for z in self._inputReader]

        # Clip the list of ZMWs if nZMWs was set
        if options.nZmws < 0:
            self._sequencingZmws = sequencingZmws
        else:
            self._sequencingZmws = sequencingZmws[:options.nZmws]

    def _loadChemistry(self):
        logging.info("Loading Chemistry information")
        all_chemistries = [self._inputReader[k].sequencingChemistry for k in self._inputReader.readers.keys()]
        used_chemistries = sorted(set(all_chemistries))
        assert len(used_chemistries) == 1
        self._chemistry = used_chemistries[0]

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
        self._barcodeNames = names
        self._barcodePairs = [(names[i], names[i+1]) for i in range(0,len(names)-1,2)]
        self._barcodePairNames = ["{0}--{1}".format(p[0], p[1]) for p in self._barcodePairs]

    @property
    def barcodeLength(self):
        return self._barcodeLength
    @property
    def windowSize(self):
        return self._windowSize
    @property
    def adapterPad(self):
        return self._adapterPad
    @property
    def insertPad(self):
        return self._insertPad
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

    #def _makeRead( self, zmw, start, end ):
    #    if (end-start) < self._barcodeLength:
    #        logging.debug("ZMW #{0} - Skipping barcode window ({1}, {2}) - too small".format(zmw.holeNumber, start, end))
    #        return None
    #    #bax = self._getBaxForHole( holeNum )
    #    return ConsensusCoreRead(zmw.baxH5, zmw.holeNumber, start, end, self.chemistry)

    def _barcodeToPairScores( self, bcScores, barcodeNames ):
        pairScores = {}
        for fwd, rev in self._barcodePairs:
            pair = "{0}--{1}".format(fwd, rev)
            fwdIdx = barcodeNames.index(fwd)
            revIdx = barcodeNames.index(rev)
            pairScores[pair] = bcScores[fwdIdx] + bcScores[revIdx]
        return pairScores

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

    def isBadAdp( self, bestIdx, errors ):
        if sum(errors) != 1:
            return False
        if errors[0] == 1 or errors[-1] == 1:
            return False
        idx = errors.index(1)
        if bestIdx[idx-1] != bestIdx[idx+1]:
            return True
        return False

    def isMissingAdp( self, bestIdx, errors):
        if len(set(bestIdx)) == 1:
            return False
        if sum(errors) != 0:
            return False
        for i in range(1, len(bestIdx)):
            if bestIdx[i] == bestIdx[i-1]:
                return True
        return False

    def isEndError( self, errors ):
        maxErrors = sum(errors)
        for i in range(1,len(errors)):
            left  = errors[:i]
            right = errors[i:]
            if len(left) == sum(left)   and sum(left) == maxErrors:
                return True
            if len(right) == sum(right) and sum(right) == maxErrors:
                return True
        return False

    def isMixedError( self, errors ):
        mixedA = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        mixedB = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        errorLen = len(errors)
        if errors == mixedA[:errorLen]:
            return True
        elif errors == mixedB[:errorLen]:
            return True
        return False

    def main(self):
        parseOptions()
        self._setupLogging()

        logging.info("h5py version: %s" % h5py.version.version)
        logging.info("hdf5 version: %s" % h5py.version.hdf5_version)
        logging.info("ConsensusCore version: %s" % consensusCoreVersion())
        logging.info("BarcodeAnalysis version: %s" % __version__)
        logging.info("Custom SO File: %s" % options.soFile)

        logging.info("Starting.")
        self._loadWhiteList()
        self._loadData()
        self._loadChemistry()
        self._loadBarcodes()

        self._insertPad = options.insertSidePad
        self._adapterPad = options.adapterSidePad

        scorer = BarcodeScorer(self._inputReader,
                               FastaReader(options.barcodeFilename),
                               adapterSidePad = options.adapterSidePad,
                               insertSidePad = options.insertSidePad,
                               scoreMode = 'paired', maxHits = options.maxHits,
                               scoreFirst = False, startTimeCutoff = 1,
                               minScore = 30,
                               soFile=options.soFile)
        # If tSNE was selected, output the barcode information as a CSV
        if options.tSNE:
            with self.openOutputFile() as handle:
                for holeNum in self._sequencingZmws:
                    zmw = self.inputReader[holeNum]
                    if len(zmw.adapterRegions) >= 3:
                        for window in self._getWindowReads(zmw):
                            handle.write( window.to_csv + '\n' )
        # If AdapterSizes was selected, output the length of the barcode on either side of the adapter
        elif options.adapterSizes:
            with self.openOutputFile() as handle:
                for zmw in self._sequencingZmws:
                    lengths = [(len(w) if w else 0) for w in self._getWindowReads(zmw)]
                    print '%s,%s,%s' % (zmw.zmwName, ','.join([str(l) for l in lengths]), min(lengths))
        # If PBbarcode was selected, use Bullard's BarcodeScorer to score sequences
        elif options.scoreBarcodes:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmw2( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxCorrect = [1 if idx in trueIdxPts else 0 for idx in adpBestIdx]
                adpScoreCorrect =   [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 1]
                adpScoreIncorrect = [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 0]
                avgCorrect   = sum(adpScoreCorrect)/float(len(adpScoreCorrect))     if len(adpScoreCorrect)   else 'N/A'
                avgIncorrect = sum(adpScoreIncorrect)/float(len(adpScoreIncorrect)) if len(adpScoreIncorrect) else 'N/A'
                print "{0},{1},{2},{3},{4},{5}".format(zmw.zmwName, trueIdx,
                                                       len(adpBestArg), sum(adpIdxCorrect),
                                                       avgCorrect, avgIncorrect)
        elif options.scoreBarcodesOld:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmw3( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxCorrect = [1 if idx in trueIdxPts else 0 for idx in adpBestIdx]
                adpScoreCorrect =   [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 1]
                adpScoreIncorrect = [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 0]
                avgCorrect   = sum(adpScoreCorrect)/float(len(adpScoreCorrect))     if len(adpScoreCorrect)   else 'N/A'
                avgIncorrect = sum(adpScoreIncorrect)/float(len(adpScoreIncorrect)) if len(adpScoreIncorrect) else 'N/A'
                print "{0},{1},{2},{3},{4},{5}".format(zmw.zmwName, trueIdx,
                                                       len(adpBestArg), sum(adpIdxCorrect),
                                                       avgCorrect, avgIncorrect)
        elif options.scoreBarcodesRc:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwRc( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxCorrect = [1 if idx in trueIdxPts else 0 for idx in adpBestIdx]
                adpScoreCorrect =   [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 1]
                adpScoreIncorrect = [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 0]
                avgCorrect   = sum(adpScoreCorrect)/float(len(adpScoreCorrect))     if len(adpScoreCorrect)   else 'N/A'
                avgIncorrect = sum(adpScoreIncorrect)/float(len(adpScoreIncorrect)) if len(adpScoreIncorrect) else 'N/A'
                print "{0},{1},{2},{3},{4},{5}".format(zmw.zmwName, trueIdx,
                                                       len(adpBestArg), sum(adpIdxCorrect),
                                                       avgCorrect, avgIncorrect)
        elif options.testBarcodesRc:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwRc( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                uniqueBestIdx = sorted(set(adpBestIdx))
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                print zmw.zmwName, adpIdxIncorrect
                if sum(adpIdxIncorrect) > 0:
                    print zmw.zmwName, trueIdx
                    print adpBestIdx
                    scorer.scoreSelectedAdaptersRc(zmw, adpIdxIncorrect, uniqueBestIdx)

        elif options.testBarcodesRc2:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                hqStart, hqEnd = zmw.hqRegion
                for adp, scores in zip(zmw.adapterRegions, scorer.scoreZmwRc2( zmw )):
                    leftEnd, rightStart = adp
                    leftScore, rightScore = scores
                    if leftEnd > hqStart:
                        leftStart = max(0, leftEnd - self.barcodeLength - self.insertPad)
                        leftSeq = reverse_complement(zmw.read(leftStart, leftEnd).basecalls())
                        leftMax = max(leftScore)
                        leftIdx = list(leftScore).index(leftMax)
                        leftBc = self._barcodeNames[leftIdx]
                        leftBcSeq = self._barcodeSequences[(leftBc, 'FORWARD')] if leftBc.startswith('F') else self._barcodeSequences[(leftBc, 'REVERSE')]
                        print "{0},{1},{2},Left,{3},{4}".format(zmw.read().readName, leftEnd, rightStart, leftBc, leftMax)
                    if rightStart < hqEnd:
                        rightEnd  = min(hqEnd, rightStart + self.barcodeLength + self.insertPad)
                        rightSeq = zmw.read(rightStart, rightEnd).basecalls()
                        rightMax = max(rightScore)
                        rightIdx = list(rightScore).index(rightMax)
                        rightBc = self._barcodeNames[rightIdx]
                        rightBcSeq = self._barcodeSequences[(rightBc, 'FORWARD')] if rightBc.startswith('F') else self._barcodeSequences[(rightBc, 'REVERSE')]
                        rightTestScore = list(rightScore)[leftIdx]
                        print "{0},{1},{2},Right,{3},{4}".format(zmw.read().readName, leftEnd, rightStart, rightBc, rightMax)

        elif options.testBarcodesRc3:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                hqStart, hqEnd = zmw.hqRegion
                for adp, scores in zip(zmw.adapterRegions, scorer.scoreZmwRc2( zmw )):
                    leftEnd, rightStart = adp
                    leftScore, rightScore = scores
                    if leftEnd > hqStart:
                        leftStart = max(0, leftEnd - self.barcodeLength - self.insertPad)
                        leftSeq = reverse_complement(zmw.read(leftStart, leftEnd).basecalls())
                        leftMax = max(leftScore)
                        leftIdx = list(leftScore).index(leftMax)
                        leftBc = self._barcodeNames[leftIdx]
                        leftBcSeq = self._barcodeSequences[(leftBc, 'FORWARD')] if leftBc.startswith('F') else self._barcodeSequences[(leftBc, 'REVERSE')]
                        print "{0},{1},{2},Left,{3},{4}".format(zmw.read().readName, leftEnd, rightStart, leftBc, leftMax)
                        scorer.aligner.score(leftSeq, leftBcSeq)
                    if rightStart < hqEnd:
                        rightEnd  = min(hqEnd, rightStart + self.barcodeLength + self.insertPad)
                        rightSeq = zmw.read(rightStart, rightEnd).basecalls()
                        rightMax = max(rightScore)
                        rightIdx = list(rightScore).index(rightMax)
                        rightBc = self._barcodeNames[rightIdx]
                        rightBcSeq = self._barcodeSequences[(rightBc, 'FORWARD')] if rightBc.startswith('F') else self._barcodeSequences[(rightBc, 'REVERSE')]
                        rightTestScore = list(rightScore)[leftIdx]
                        print "{0},{1},{2},Right,{3},{4}".format(zmw.read().readName, leftEnd, rightStart, rightBc, rightMax)
                        scorer.aligner.score(rightSeq, leftBcSeq)

        elif options.funnyAdapters:
            for zmw in self._sequencingZmws:
                adpEnds = [str(r[1]) for r in zmw.adapterRegions[:options.maxHits]]
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwRc( zmw, trim=True )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                isBadAdp = self.isBadAdp( adpBestIdx, adpIdxIncorrect )
                isMissingAdp = self.isMissingAdp( adpBestIdx, adpIdxIncorrect )
                bestIdxStr = '--'.join(adpBestIdx)
                adpEndStr = '--'.join(adpEnds)
                if isBadAdp:
                    print "{0},ExtraAdp,{1},{2}".format(zmw.zmwName, bestIdxStr, adpEndStr)
                elif isMissingAdp:
                    print "{0},MissingAdp,{1},{2}".format(zmw.zmwName, bestIdxStr, adpEndStr)
                elif sum(adpIdxIncorrect) == 0 and not isMissingAdp:
                    print "{0},Perfect,{1},{2}".format(zmw.zmwName, bestIdxStr, adpEndStr)

        elif options.summarizeErrorsRc:
            errors = {'total': 0, 'true': 0, 'mixed':0, 'other': 0, 'badAdp':0, 'single':0, 'multi':0, 'map':0, 'bad':0, 'readScore':0}
            for zmw in self._sequencingZmws:
                # Score the Zmw as per normal
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwRc( zmw )
                adpScores = res[1]

                # Figure out the best scores and Ids for each Adapter
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]

                # Figure out which were correct, and which are in error
                adpIdxCorrect = [1 if idx in trueIdxPts else 0 for idx in adpBestIdx]
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                numErrors = sum(adpIdxIncorrect)

                # Test for certain error patterns
                isBadAdp = self.isBadAdp( adpBestIdx, adpIdxIncorrect )
                isEndError = self.isEndError( adpIdxIncorrect )
                isMixedError = self.isMixedError( adpIdxIncorrect )

                # Figure out the scores and averages for the correct/incorrect barcode calls
                adpScoreCorrect =   [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 1]
                adpScoreIncorrect = [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 0]
                avgCorrect   = sum(adpScoreCorrect)/float(len(adpScoreCorrect))     if len(adpScoreCorrect)   else 'N/A'
                avgIncorrect = sum(adpScoreIncorrect)/float(len(adpScoreIncorrect)) if len(adpScoreIncorrect) else 'N/A'

                # Summarize the errors
                errors['total'] += numErrors
                if numErrors == 1:
                    errors['single'] += numErrors
                else:
                    errors['multi'] += numErrors

                if avgCorrect == 'N/A':
                    errors['map'] += numErrors
                elif isBadAdp:
                    errors['badAdp'] += numErrors
                elif isEndError and avgIncorrect > 30.0:
                    errors['true'] += numErrors
                elif isMixedError:
                    errors['mixed'] += numErrors
                elif numErrors == 1 and avgIncorrect < 30.0:
                    errors['bad'] += numErrors
                elif zmw.readScore <= 0.8:
                    errors['readScore'] += numErrors
                else:
                    errors['other'] += numErrors
                    print zmw.zmwName
                    print adpBestIdx
                    print adpIdxIncorrect
                    print avgCorrect, avgIncorrect
                    print isBadAdp, isEndError, isMixedError

            print errors

        # Otherwise score the barcodes normally and return a CSV
        elif options.pbbarcode2:
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmw( zmw )
                adpScores = res[1]
                print scorer.barcodeNames
                print res[0]
                print res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                print adpBestArg
                print [max(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                uniqueBestIdx = sorted(set(adpBestIdx))
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                print zmw.zmwName
                print adpBestIdx
                print adpIdxIncorrect
                scorer.scoreSelectedAdapters(zmw, adpIdxIncorrect, uniqueBestIdx)
                print
        elif options.pbbarcode3:
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmw2( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpTrueScores = [v for i,v in enumerate(adpBestScores) if adpIdxIncorrect[i]==0]
                adpTrueScores += (10-len(adpTrueScores)) * [0]
                scoresAsStr = [str(s) for s in adpTrueScores]
                print "{0},{1},{2}".format(zmw.zmwName, trueIdx, ",".join(scoresAsStr))
                #try:
                #    firstCorrect = adpIdxIncorrect.index(0)
                #    scores = sorted(adpScores[firstCorrect])
                #    scoresAsStr = [str(s) for s in scores]
                #    print "{0},{1},{2},{3}".format(zmw.zmwName, trueIdx, ",".join(scoresAsStr[:-1]))
                #except:
                #    continue
        # Otherwise score the barcodes normally and return a CSV
        elif options.scoreAdapters:
            print "Zmw,TrueIdx,NumAdp,NumCorrect,CorrectAvg,IncorrectAvg"
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwAdps( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                adpIdxCorrect = [1 if idx in trueIdxPts else 0 for idx in adpBestIdx]
                adpScoreCorrect =   [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 1]
                adpScoreIncorrect = [adpBestScores[i] for i,v in enumerate(adpIdxCorrect) if v == 0]
                avgCorrect   = sum(adpScoreCorrect)/float(len(adpScoreCorrect))     if len(adpScoreCorrect)   else 'N/A'
                avgIncorrect = sum(adpScoreIncorrect)/float(len(adpScoreIncorrect)) if len(adpScoreIncorrect) else 'N/A'
                print "{0},{1},{2},{3},{4},{5}".format(zmw.zmwName, trueIdx,
                                                       len(adpBestArg), sum(adpIdxCorrect),
                                                       avgCorrect, avgIncorrect)
        elif options.testAdapters:
            for zmw in self._sequencingZmws:
                trueIdx = self._whiteListIdx[ zmw.zmwName ]
                trueIdxPts = trueIdx.split('--')
                res = scorer.scoreZmwAdps( zmw )
                adpScores = res[1]
                adpBestArg = [np.argmax(a) for a in adpScores]
                adpBestIdx = [scorer.barcodeNames[i] for i in adpBestArg]
                uniqueBestIdx = sorted(set(adpBestIdx))
                adpIdxIncorrect = [0 if idx in trueIdxPts else 1 for idx in adpBestIdx]
                adpBestScores = [adpScores[i][x] for i,x in enumerate(adpBestArg)]
                if sum(adpIdxIncorrect) > 0:
                    print zmw.zmwName, trueIdx
                    print adpBestIdx
                    scorer.scoreSelectedAdapterRegions(zmw, adpIdxIncorrect, uniqueBestIdx)
        else:
            barcodeNames = scorer.barcodeNames2
            with self.openOutputFile() as handle:
                for zmw in self._sequencingZmws:
                    if len(zmw.adapterRegions) < 1:
                        continue
                    trueIdx = self._whiteListIdx[ zmw.zmwName ]
                    res = scorer.scoreZmwRc( zmw )
                    bcScores = res[0]
                    pairScores = self._barcodeToPairScores( bcScores, barcodeNames )
                    if bcScores is not None:
                        handle.write( barcodeCsvLine2( zmw, trueIdx, pairScores ) )

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
