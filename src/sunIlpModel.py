import pulp
from collections import OrderedDict, defaultdict
import numpy as np
from itertools import izip

from src.abstractIlpSolving import SequenceGraphLpProblem, get_id
from lib.general_lib import rejectOutliers


class Window(object):
    """Stores the allele fractions of the A and B probes in a specific window
    as well as the LP variables to be solved for this window"""
    def __init__(self, A, B):
        A_values = rejectOutliers(np.array(A))
        B_values = rejectOutliers(np.array(B))
        if len(A_values) == 0:
            self.A_value = 2
        else:
            self.A_value = np.mean(A_values)
        if len(B_values) == 0:
            self.B_value = 2
        else:
            self.B_value = np.mean(B_values)
        #save the LP variables that represent the inferred copy number at this window
        self.A = pulp.LpVariable("A_{}".format(get_id()), lowBound=0, upBound=4, cat="Integer")
        self.B = pulp.LpVariable("B_{}".format(get_id()), lowBound=0, upBound=4, cat="Integer")


class SunIlpModel(SequenceGraphLpProblem):
    def __init__(self, Avals, Bvals, windowSize, stepSize, breakpoint_penalty, data_penalty, deletion_penalty):
        aStart = 146152644-3000
        bStart = 148603586-3000
        aStop = 146233816+4000
        bStop = 148684557+4000
        SequenceGraphLpProblem.__init__(self)
        self.windows = []
        for aPos, bPos in izip(xrange(aStart, aStop - windowSize, stepSize), xrange(bStart, bStop - windowSize, stepSize)):
            A = [y for x, y in Avals if x >= aPos and y < aPos]
            B = [y for x, y in Bvals if x >= bPos and y < bPos]
            self.windows.append([aPos, bPos, Window(A, B)])
        self.build_model(breakpoint_penalty, data_penalty, deletion_penalty)

    def get_results(self):
        """generator that yields mid, Aval, Bval results. ILP must be solved"""
        assert self.is_solved
        for aPos, bPos, window in self.windows:
            yield aPos, bPos, pulp.value(window.A), pulp.value(window.B)

    def build_model(self, breakpoint_penalty, data_penalty, deletion_penalty):
        """Builds ILP model. Run this once before solving the model.
        See the program docstring for the math behind this"""
        for i in xrange(len(self.windows)):
            window = self.windows[i][-1]
            #minimize differences between data and variables
            self.constrain_approximately_equal(window.A_value, window.A, data_penalty)
            self.constrain_approximately_equal(window.B_value, window.B, data_penalty)
            #minimize deviation from 2 copies of A and B
            self.constrain_approximately_equal(sum([window.A, window.B]), 4, deletion_penalty)
            if i != len(self.windows) - 1:
                #penalize introducing breakpoints; tie windows together
                next_window = self.windows[i + 1][-1]
                self.constrain_approximately_equal(window.A, next_window.A, breakpoint_penalty)
                self.constrain_approximately_equal(window.B, next_window.B, breakpoint_penalty)

    def run(self):
        self.solve()
        results = []
        for aPos, bPos, aResult, bResult in self.get_results():
            results.append([aPos, aResult])
            results.append([bPos, bResult])
        return sorted(results, key = lambda x: x[0])