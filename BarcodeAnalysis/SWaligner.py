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
from ctypes import *
import os
import numpy
import pkg_resources

class SWaligner(object):
    def __init__(self, soFile=None):
        # setup.py should put sw.so in the following path.
        if soFile is None:
            self.SW_DLL_PATH = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "sw.so"
        else:
            self.SW_DLL_PATH = os.path.abspath( soFile )
        self._dll        = CDLL(self.SW_DLL_PATH)
        self.dpMat       = self._dll.allocate_dp_mat()

    def score(self, tSeq, qSeq):
        score = self._dll.compute_align_score(self.dpMat, tSeq, qSeq)
        self._dll.print_dp_mat(self.dpMat, tSeq, qSeq)
        print "Max: %s" % score
        return score

    def makeScorer(self, targets):
        ScoreType = c_int * len(targets)
        scores = ScoreType()
        for i in range(0, len(scores)):
            scores[i] = 0

        TargetType = c_char_p * len(targets)
        targetSeqs = TargetType()
        for i in range(0, len(targetSeqs)):
            targetSeqs[i] = targets[i]

        targetLen = len(targets)

        def scorer(query):
            if not query:
                return numpy.zeros(len(targets))

            self._dll.compute_align_scores(scores,
                                           targetLen,
                                           self.dpMat,
                                           query,
                                           targetSeqs)
            return numpy.array([scores[i] for i in xrange(0, len(scores))])

        return scorer
