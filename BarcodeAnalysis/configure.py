from __future__ import absolute_import

import ConfigParser

import ConsensusCore as cc

from BarcodeAnalysis.utils import fst, snd

__author__ = 'Brett Bowman'

BANDING_DIAG = 4
BANDING_SCORE = 50
FAST_SCORE_THRESH = -12.5

def _buildQuiverConfig(parameterSetName, nameValuePairs):
    qvModelParams      = cc.QvModelParams(*[ float(snd(pair)) for pair in nameValuePairs ])
    bandingOptions     = cc.BandingOptions( BANDING_DIAG, BANDING_SCORE )
    fastScoreThreshold = FAST_SCORE_THRESH
    return  cc.QuiverConfig(qvModelParams,
                            cc.ALL_MOVES,
                            bandingOptions,
                            fastScoreThreshold)

def loadQuiverConfigs(iniFilename):
    # returns dict: name -> ParameterSet
    cp = ConfigParser.ConfigParser()
    cp.optionxform=str
    cp.read([iniFilename])
    sections = cp.sections()
    quiverConfigs = {}
    for sectionName in sections:
        config = _buildQuiverConfig(sectionName, cp.items(sectionName))
        if config:
            quiverConfigs[sectionName] = config
    return quiverConfigs
