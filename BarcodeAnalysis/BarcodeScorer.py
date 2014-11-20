#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$
import logging

from pbcore.io import BasH5Reader, BaxH5Reader
from pbcore.io.FastaIO import *
import BarcodeAnalysis.SWaligner as Aligner
import numpy as n

from pbcore.io.BarcodeH5Reader import LabeledZmw, \
    BARCODE_DELIMITER

__RC_MAP__ = dict(zip('ACGTacgt-N','TGCAtgca-N'))

class BarcodeScorer(object):
    def __init__(self, basH5, barcodeFasta,
                 adapterSidePad = 4, insertSidePad = 4,
                 scoreMode = 'paired', maxHits = 10,
                 scoreFirst = False, startTimeCutoff = 1,
                 minScore = 20,
                 soFile=None):
        """A BarcodeScorer object scores ZMWs and produces summaries
        of the scores. Various parameters control the behavior of the
        object, specifically the padding allows the user to add a
        little extra on each side of the adapter find for safety. The
        most relevant parameter is the scoreMode which dictates how
        the barcodes are scored, either paired or symmetric."""

        self.basH5 = basH5
        self.barcodeFasta = list(barcodeFasta)
        self.aligner = Aligner.SWaligner(soFile)
        self.barcodeLength = n.unique(map(lambda x : len(x.sequence),
                                          self.barcodeFasta))
        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

        self.adapterSeq   = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"
        self.adapterSeqRc = self._rc("ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT")
        self.barcodeSeqs = [(barcode.sequence.upper(),
                             self._rc(barcode.sequence.upper()))
                            for barcode in self.barcodeFasta]

        self.adapterSidePad = adapterSidePad
        self.insertSidePad = insertSidePad
        self.maxHits = maxHits
        self.minScore = minScore

        if scoreMode not in ['symmetric', 'paired']:
            raise Exception("scoreMode must either be symmetric or paired")
        self._scoreMode = scoreMode

        self.scoreFirst = scoreFirst
        self.startTimeCutoff = startTimeCutoff

        self.threePrimeSeqs  = [(x[0] if (i%2)==0 else x[1]) for i,x in enumerate(self.barcodeSeqs)]
        self.fivePrimeSeqs   = [(x[1] if (i%2)==0 else x[0]) for i,x in enumerate(self.barcodeSeqs)]
        self.fivePrimeSeqsRc = [self._rc(s) for s in self.fivePrimeSeqs]
        def fwdBcAdp( seq_pair ):
            return "{0}{1}{2}".format(seq_pair[1], self.adapterSeqRc, seq_pair[0])
        def revBcAdp( seq_pair ):
            return "{0}{1}{2}".format(seq_pair[0], self.adapterSeqRc, seq_pair[1])
        self.adapterBcSeqs = [fwdBcAdp(p) if (i%2)==0 else revBcAdp(p) for i,p in enumerate(self.barcodeSeqs)]

        self.threePrimeScorer  = self.aligner.makeScorer(self.threePrimeSeqs)
        self.fivePrimeScorer   = self.aligner.makeScorer(self.fivePrimeSeqs)
        self.fivePrimeScorerRc = self.aligner.makeScorer(self.fivePrimeSeqsRc)
        self.adapterScorer     = self.aligner.makeScorer(self.adapterBcSeqs)

        self.forwardScorer = self.aligner.makeScorer([x[0] for x in self.barcodeSeqs])
        self.reverseScorer = self.aligner.makeScorer([x[1] for x in self.barcodeSeqs])

        logging.debug(("Constructed BarcodeScorer with scoreMode: %s," + \
                           "adapterSidePad: %d, insertSidePad: %d, and scoreFirst: %r") \
                          % (scoreMode, adapterSidePad, insertSidePad, scoreFirst))

    @property
    def movieName(self):
        return self.basH5.movieName

    def makeBCLabel(self, s1, s2):
        return BARCODE_DELIMITER.join((s1, s2))

    @property
    def barcodeLabels(self):
        """The barcode labels are function of the barcodeNames and the
        scoreMode, they represent the user-visible names."""
        if self.scoreMode == 'paired':
            return n.array([self.makeBCLabel(self.barcodeFasta[i].name,
                                             self.barcodeFasta[i+1].name) for i
                            in xrange(0, len(self.barcodeSeqs), 2)])
        else:
            return n.array([self.makeBCLabel(x.name, x.name) for x in self.barcodeFasta])

    @property
    def barcodeNames(self):
        """The barcode names are the FASTA names"""
        return n.array([x.name for x in self.barcodeFasta])

    @property
    def barcodeNames2(self):
        return [x.name for x in self.barcodeFasta]

    @property
    def scoreMode(self):
        return self._scoreMode

    def _rc(self, s):
        return "".join([__RC_MAP__[c] for c in s[::-1]])

    def _adapterSeqs(self, zmw):
        def fromRange(rStart, rEnd):
            try:
                adpSeq = zmw.read(rStart - (self.barcodeLength + self.insertSidePad),
                                  rEnd   +  self.barcodeLength + self.insertSidePad).basecalls()
            except IndexError:
                return None
            return adpSeq

        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]

        seqs = [fromRange(start, end) for (start, end) in adapterRegions]
        return seqs

    def _flankingSeqs(self, zmw):
        def fromRange(rStart, rEnd):
            try:
                qSeqLeft = zmw.read(rStart - (self.barcodeLength + self.insertSidePad),
                                    rStart + self.adapterSidePad).basecalls()
            except IndexError:
                qSeqLeft = None
            try:
                qSeqRight = zmw.read(rEnd - self.adapterSidePad,
                                     rEnd + self.barcodeLength +
                                     self.insertSidePad).basecalls()
            except IndexError:
                qSeqRight = None

            return (qSeqLeft, qSeqRight)

        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]

        seqs = [fromRange(start, end) for (start, end) in adapterRegions]

        # We only score the first barcode if we don't find any adapters
        # *and* the start time is less than the threshold.
        scoredFirst = False
        if self.scoreFirst and not len(seqs):
            s = zmw.zmwMetric('HQRegionStartTime')
            e = zmw.zmwMetric('HQRegionEndTime')
            # s<e => has HQ.
            if s < e and s <= self.startTimeCutoff:
                l = self.barcodeLength + self.insertSidePad
                l = l if zmw.hqRegion[1] > l else zmw.hqRegion[1]
                try:
                    bc = zmw.read(0, l).basecalls()
                    if len(bc) >= self.barcodeLength:
                        seqs.insert(0, (bc, None))
                        scoredFirst = True
                except IndexError:
                    pass

        return (seqs, scoredFirst)

    def testAligner(self, holeNumbers):
        for holeNumber in holeNumbers:
            print holeNumber
            zmw = self.basH5[holeNumber]
            print zmw
            adapters, _ = self._flankingSeqs(zmw)
            for left, right in adapters:
                for fwd, rev in self.barcodeSeqs:
                    print len(fwd), len(rev)
                    if left:
                        print "Left, Fwd"
                        self.aligner.score(left, fwd)
                        print "Left, Rev"
                        self.aligner.score(left, rev)
                    if right:
                        print "Right, Fwd"
                        self.aligner.score(right, fwd)
                        print "Right, Rev"
                        self.aligner.score(right, rev)

    def scoreZmw(self, zmw):
        adapters, scoredFirst = self._flankingSeqs(zmw)
        adapterScores = [[]]*len(adapters)
        barcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            fscores  = self.forwardScorer(adapter[0])
            rscores  = self.reverseScorer(adapter[0])
            ffscores = self.forwardScorer(adapter[1])
            rrscores = self.reverseScorer(adapter[1])

            scored = 2.0 if adapter[0] and adapter[1] \
                else 1.0 if adapter[0] or  adapter[1] \
                else 0

            # An adapter score is the average barcode score for
            # each barcode -- that way, you can compare across
            # adapters even if the different adapters have
            # different numbers of flanking sequence.
            if scored == 0:
                adapterScores[i] = barcodeScores
            else:
                adapterScores[i] = n.maximum((fscores + rrscores)/scored,
                                             (rscores + ffscores)/scored)

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else n.zeros(len(self.barcodeSeqs))

        return (barcodeScores, adapterScores)

    def scoreZmw2(self, zmw):
        adapters, scoredFirst = self._flankingSeqs(zmw)
        adapterScores = [[]]*len(adapters)
        barcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            fscores = self.fivePrimeScorer(adapter[0])
            tscores = self.threePrimeScorer(adapter[1])

            filteredF      = n.array([(s if s >= self.minScore else 0) for s in fscores])
            filteredT      = n.array([(s if s >= self.minScore else 0) for s in tscores])
            #filteredCounts = n.array([(2.0 if filteredF[i] and filteredT[i] else 1.0) for i in range(len(fscores))])

            scored = 2.0 if adapter[0] and adapter[1] else \
                1.0 if adapter[0] or adapter[1] else 0

            # An adapter score is the average barcode score for
            # each barcode -- that way, you can compare across
            # adapters even if the different adapters have
            # different numbers of flanking sequence.
            if scored == 0:
                adapterScores[i] = barcodeScores
            else:
                adapterScores[i] = (filteredF + filteredT)/scored
                #adapterScores[i] = (fscores + tscores)/scored

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else n.zeros(len(self.barcodeSeqs))

        return (barcodeScores, adapterScores)

    def scoreZmw3(self, zmw):
        adapters, scoredFirst = self._flankingSeqs(zmw)
        adapters2 = [((a[0], a[1]) if a[0] is None else (self._rc(a[0]), a[1])) for a in adapters]
        adapterScores = [[]]*len(adapters)
        barcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            fscores = self.fivePrimeScorer2(adapter[0])
            tscores = self.threePrimeScorer(adapter[1])

            scored = 2.0 if adapter[0] and adapter[1] \
                else 1.0 if adapter[0] or adapter[1]  \
                else 0

            # An adapter score is the average barcode score for
            # each barcode -- that way, you can compare across
            # adapters even if the different adapters have
            # different numbers of flanking sequence.
            if scored == 0:
                adapterScores[i] = barcodeScores
            else:
                #adapterScores[i] = (filteredF + filteredT)/scored
                adapterScores[i] = (fscores + tscores)/scored

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else n.zeros(len(self.barcodeSeqs))

        return (barcodeScores, adapterScores)

    def scoreZmwAdps(self, zmw):
        adapters = self._adapterSeqs(zmw)
        perAdapterScores = [[]]*len(adapters)
        perBarcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            perAdapterScores[i] = self.adapterScorer(adapter)

        perBarcodeScores = reduce(lambda x, y: x + y, perAdapterScores) \
                                      if perAdapterScores \
                                      else n.zeros(len(self.barcodeSeqs))

        return (perBarcodeScores, perAdapterScores)

    def scoreSelectedAdapters(self, zmw, selectedAdp, selectedBc):
        adapters, scoredFirst = self._flankingSeqs(zmw)
        assert len(adapters) == len(selectedAdp)
        selectedAdapters = [adapters[i] for i,v in enumerate(selectedAdp) if v == 1]
        selectedBcPos = [self.barcodeNames2.index(bc) for bc in selectedBc]
        selectedBcSeqs = [self.barcodeSeqs[i] for i in selectedBcPos]
        for i, adps in enumerate(selectedAdapters):
            fwdAdp, revAdp = adps
            print "FORWARD"
            for j, bc in enumerate(selectedBc):
                fwdBc, revBc = selectedBcSeqs[j]
                print "Adp #{0} - BC {1} - FwdAdp FwdBc".format(i+1, bc)
                self.aligner.score(fwdAdp, fwdBc)
                print "Adp #{0} - BC {1} - FwdAdp RevBc".format(i+1, bc)
                self.aligner.score(fwdAdp, revBc)
            print "REVERSE"
            for j, bc in enumerate(selectedBc):
                fwdBc, revBc = selectedBcSeqs[j]
                print "Adp #{0} - BC {1} - RevAdp FwdBc".format(i+1, bc)
                self.aligner.score(revAdp, fwdBc)
                print "Adp #{0} - BC {1} - RevAdp RevBc".format(i+1, bc)
                self.aligner.score(revAdp, revBc)
            print "END\n"

    def scoreSelectedAdapterRegions(self, zmw, selectedAdp, selectedBc):
        adapters = self._adapterSeqs(zmw)
        assert len(adapters) == len(selectedAdp)
        selectedAdapters = [adapters[i] for i,v in enumerate(selectedAdp) if v == 1]
        selectedBcPos = [self.barcodeNames2.index(bc) for bc in selectedBc]
        selectedAdpBcSeqs = [self.adapterBcSeqs[i] for i in selectedBcPos]
        print zmw.zmwName
        for i, adp in enumerate(selectedAdapters):
            for j, bcId in enumerate(selectedBc):
                adpBc = selectedAdpBcSeqs[j]
                print "Adp #{0} - BC {1} - FwdAdp FwdBc".format(i+1, bcId)
                self.aligner.score(adp, adpBc)
            print "END\n"

    def chooseSymmetric(self, o):
        p = n.argsort(-o[2])
        return LabeledZmw(o[0], o[1], p[0], o[2][p[0]], p[1], o[2][p[1]], o[3])

    def choosePaired(self, o):
        if o[1] == 1:
            s = n.array([max(o[2][i], o[2][i + 1]) for i in \
                             xrange(0, len(self.barcodeSeqs), 2)])
            p = n.argsort(-s)
            s = s[p]
        else:
            # score the pairs by scoring the two alternate
            # ways they could have been put on the molecule. A
            # missed adapter will confuse this computation.
            scores  = o[3]
            results = n.zeros(len(self.barcodeSeqs)/2)
            for i in xrange(0, len(self.barcodeSeqs), 2):
                pths = [0,0]
                for j in xrange(0, len(scores)):
                    pths[j % 2] += scores[j][i]
                    pths[1 - j % 2] += scores[j][i + 1]
                results[i/2] = max(pths)

            p = n.argsort(-results)
            s = results[p]

        return LabeledZmw(o[0], o[1], p[0], s[0], p[1], s[1], o[3])

    def labelZmws(self, holeNumbers):
        """Return a list of LabeledZmws for input holeNumbers"""
        # o here is the record immediately above.
        if self.scoreMode == 'symmetric':
            choose = self.chooseSymmetric
        elif self.scoreMode == 'paired':
            choose = self.choosePaired
        else:
            raise Exception("Unsupported scoring mode in BarcodeLabeler.py")

        scored = [self.scoreZmw(self.basH5[zmw]) for zmw in holeNumbers]
        return [choose(scoreTup) for scoreTup in scored if scoreTup[1]]
