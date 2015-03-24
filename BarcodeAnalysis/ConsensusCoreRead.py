__author__ = "Brett Bowman"
"""
import ConsensusCore as cc

from BarcodeAnalysis.utils import (arrayFromDataset,
                                   asFloatFeature,
                                   QUIVER_FEATURES)

class ConsensusCoreRead(object):
    def __init__(self, bax, holeNum, start, end, chemistry):
        self._bax = bax
        self._holeNum = holeNum
        self._offset = self._getOffset()
        self._start = start
        self._end = end
        self._chemistry = chemistry
        self._sequence = self._getSequence()
        self._qvs = self._getQvFeatures()
        self._sequenceFeatures = self._getSequenceFeatures()
        self._read = cc.Read(self._sequenceFeatures, self.name, self.chemistry)

    def _getOffset(self):
        return self._bax._offsetsByHole[self._holeNum][0]

    def _getSequenceFeature( self, feature ):
        return arrayFromDataset(self._bax._basecallsGroup[feature],
                                self.absStart, self.absEnd)

    def _getQvFeatures(self):
        qvs = {}
        for feature in QUIVER_FEATURES:
            qvs[feature] = self._getSequenceFeature( feature )
        return qvs

    def _getSequence(self):
        return self._getSequenceFeature( "Basecall" ).tostring()

    def _getSequenceFeatures(self):
        features = [self.sequence]
        for feature in QUIVER_FEATURES:
            features.append( asFloatFeature( self._qvs[feature] ) )
        return cc.QvSequenceFeatures(*features)

    @property
    def movieName(self):
        return self._bax.movieName
    @property
    def absStart(self):
        return self._offset + self._start
    @property
    def absEnd(self):
        return self._offset + self._end
    @property
    def name(self):
        return "{0}/{1}/{2}_{3}".format(self.movieName, self._holeNum, self._start, self._end)
    @property
    def sequence(self):
        return self._sequence
    @property
    def chemistry(self):
        return self._chemistry
    @property
    def read(self):
        return self._read

    def _to_csv_segment(self, pos):
        segment = ',' + self._sequence[pos]
        for feature in QUIVER_FEATURES:
            segment += ',' + str(self._qvs[feature][pos])
        return segment

    def __len__(self):
        return len(self.sequence)

    @property
    def to_csv(self):
        line = self.name
        for i in range(len(self._sequence)):
            line += self._to_csv_segment(i)
        return line

    def _makeSequenceFeatures( self, holeNum, start, end ):
        absStart = zmwOffsetStart + start
        absEnd = zmwOffsetStart + end
        # Initialize the Seq features with the raw sequence
        sequenceFeatures = [arrayFromDataset(bax._basecallsGroup["Basecall"],
                            absStart, absEnd).tostring()]
        # Add each feature from the required feature set
        for feature in REQ_FEATURES:
            arrayFeature = arrayFromDataset(bax._basecallsGroup[feature], absStart, absEnd)
            floatFeature = asFloatFeature( arrayFeature )
            sequenceFeatures.append( floatFeature )
"""
