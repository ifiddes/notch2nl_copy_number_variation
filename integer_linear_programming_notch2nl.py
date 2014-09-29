import argparse, logging, sys, os, collections, math, random, itertools
import warnings, traceback, multiprocessing
import vcf
import pulp
import scipy.stats, scipy.optimize
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#hard coded region we are looking at (hg19 reference)
region=(120554845,120572645)
#new region
#region=(120554845,120579535)

class PenaltyTree(object):
    """
    Maintains a tree of penalty terms, so that we can have arbitrarily many
    penalty terms without arbitrarily large constraint expressions.
    
    """
    
    def __init__(self, degree=100):
        """
        Make a new PenaltyTree. degree specifies the maximum number of terms to
        sum together at once. Only one PenaltyTree may be used on a given LP
        problem.
        
        """
        # This holds the number of children per node/number of terms to sum at
        # once.
        self.degree = degree
        # This holds all our leaf-level terms.
        self.terms = []
        
    def get_variable(self):
        """
        Return a fresh LpVariable with a unique name.
        
        """
        # Make the variable
        var = pulp.LpVariable("PenaltyTree_{}".format(get_id()))
        # Give the fresh variable to the caller.
        return var
        
    def add_term(self, term):
        """
        Add the given LP expression as a term in the tree.
        
        """
        self.terms.append(term)
        
    def set_objective(self, problem):
        """
        Add the sum of all terms as the given LP problem's objective. The
        PenaltyTree must have at least one term.
        
        """
        # Algorithm: Go through our leaves, making a variable for the sum of
        # each group of self.degree terms. Then repeat on the sums, making sums
        # of sums, and so on. When we only have one sum, make that the
        # objective.
        # This holds the list we're collecting
        collecting = self.terms
        # This holds the list of sum variables
        sums = []
        while len(collecting) > 1:
            logging.debug("Collecting {} terms in groups of {}".format(len(
                collecting), self.degree))
            for i in xrange(0, len(collecting), self.degree):
                # This holds the terms we have collected for this sum
                collected = []
                for j in xrange(0, min(self.degree, len(collecting) - i)):
                    # Grab each term up to our degree that actually exists
                    collected.append(collecting[i + j])
                # This holds the variable we use to represent the sum of the
                # things we have collected
                sum_var = self.get_variable()
                # Constrain this variable to equal the sum of all the collected
                # terms
                problem += sum(collected) == sum_var
                # Add this variable for the next level of collection
                sums.append(sum_var)
            # Move up a level in the tree
            collecting = sums
            sums = []
        # We have now collected everything down to one term, which is in the
        # collecting list. Use it as our objective function.
        problem += collecting[0]
        
class SequenceGraphLpProblem(object):
    """
    Represents an LP copy number problem. You can attach several models to them,
    constrain them together, and solve the problem.
    
    Internally, contains a pulp LpProblem, and a PenaltyTree.
    
    """
    def __init__(self):
        """
        Make a new SequenceGraphLpProblem that we can solve.
        
        """
        # We need an actual LpProblem
        self.problem = pulp.LpProblem("copynumber", pulp.LpMinimize)
        # We also need a PenaltyTree for organizing penalty terms
        self.penalties = PenaltyTree()

    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints to the CopyNumberLpProblem's internal
        LpProblem, and penalties to the model's PenaltyTree.
        
        """
        # Make an LP variable for the amount that var_b is above var_a. Note
        # that str() on variables produces their names. Also, we have to make
        # sure that this starts with a letter.
        amount_over = pulp.LpVariable("over_{}".format(get_id()), 0)
        # Add the constraint for not being more than that much over
        self.add_constraint(var_b <= var_a + amount_over)
        # Make an LP variable for the amount that var_b is below var_a
        amount_under = pulp.LpVariable("under_{}".format(get_id()), 0)
        # Add the constraint for not being more than that much under
        self.add_constraint(var_b >= var_a - amount_under)
        # Apply an equal penalty in each direction
        self.add_penalty((penalty * amount_over) + (penalty * amount_under)) 

    def add_penalty(self, term):
        """
        Add the given penalty term to the problem's objective.
        
        """
        # Just put the term in the PenaltyTree
        self.penalties.add_term(term)

    def add_constraint(self, constraint):
        """
        Add the given (exact) constraint to the problem. For approximate
        constraints, use constrain_approximately_equal() instead.
        
        """
        # Just add the constraint to the internal problem.
        self.problem += constraint

    def solve(self, save=None):
        """
        Solve the LP problem with GLPK
        
        If save is specified, it is a filename to which to save the LP problem
        in LP format.
        
        You may only solve a SequenceGraphLpProblem once.
        
        """
        # Set up the penalties described by the penalty tree
        logging.info("Setting up penalties")
        self.penalties.set_objective(self.problem)
        if save is not None:
            logging.info("Saving problem to {}".format(save))
            self.problem.writeLP(save)
        # Solve the problem
        status = self.problem.solve(pulp.GLPK())
        logging.info("Solution status: {}".format(pulp.LpStatus[status]))        
        if len(self.problem.variables()) < 20:
            # It's short enough to look at.
            for var in self.problem.variables():
                logging.debug("\t{} = {}".format(var.name, pulp.value(var)))
        # Report the total penalty we got when we solved.
        logging.info("Penalty: {}".format(pulp.value(self.problem.objective)))
        if status != pulp.constants.LpStatusOptimal:
            raise Exception("Unable to solve problem optimally.")


class Window(object):
    """Stores the allele fractions of the A, B and C probes in a specific window
    as well as the LP variables to be solved for this window"""
    def __init__(self, A, B, pos, min_ploidy=0):
        #save the pospoint of the window
        self.pos = pos
        #save the lists of values as numpy arrays
        self.A_values = np.array(A)
        self.B_values = np.array(B)
        #save the LP variables that represent the inferred copy number at this window
        self.A = pulp.LpVariable("A_{}".format(pos), min_ploidy, cat="Integer")
        self.B = pulp.LpVariable("B_{}".format(pos), min_ploidy, cat="Integer")

    def get_values(self):
        """returns the value lists without removing outliers"""
        return [self.A_values, self.B_values]

    def get_best_estimate(self):
        """Returns a single value for A, B and C in this window that is a best estimate of A+B+C"""
        A = np.mean(reject_outliers(self.A_values))
        B = np.mean(reject_outliers(self.B_values))
        return A, B

    def get_values_outliers_removed(self):
        """Returns values with outliers removed"""
        return [reject_outliers(self.A_values), reject_outliers(self.B_values)]

    def get_copy_number(self):
        """returns copy numbers of A, B and C. LP must be solved"""
        return [pulp.value(self.A), pulp.value(self.B)]


class Model(object):
    """
    Represents model of the copy number of two genes.
    Contains a list of Window objects. Approximate-equality constraints
    can be generated tying Window copy numbers to those of their neighbors.
    
    All constraints and penalty terms are automatically added to the
    SequenceGraphLpProblem to which the Model is attached (specified on
    construction).
    
    """
    def __init__(self, problem):
        self.problem = problem
        self.windows = dict()

    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints and penalties to the Model's
        SequenceGraphLpProblem.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        # Just forward the call on to the problem.
        self.problem.constrain_approximately_equal(var_a, var_b, 
            penalty=penalty)    

    def add_constraint(self, constraint):
        """
        Add the given constraint to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        self.problem.add_constraint(constraint)

    def add_penalty(self, penalty): 
        """
        Add the given penalty term to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        self.problem.add_penalty(penalty)

    def add_window(self, window):
        """Adds a window object to the windows dict, mapped by pospoint"""
        self.windows[window.pos] = window

    def get_windows(self):
        """generator that yields pairs of pospoint, window object"""
        for pos, window in self.windows.iteritems():
            yield pos, window

    def sort_windows(self):
        """Sorts windows by pospoint"""
        self.windows = collections.OrderedDict(sorted(self.windows.items()))

    def build_model(self, default_ploidy=2, breakpoint_penalty=2, data_penalty=5, end_penalty=1):
        """Builds ILP model. Run this once before solving the model.
        See the program docstring for the math behind this"""
        self.sort_windows()
        windows = self.windows.values()
        for i in xrange(len(windows)):
            window = windows[i]
            #find best data values in this window
            D_A, D_B = window.get_best_estimate()

            #minimize differences between data and variables
            self.constrain_approximately_equal(D_A, window.A, data_penalty)
            self.constrain_approximately_equal(D_B, window.B, data_penalty)

            #penalize introducing breakpoints; tie windows together
            if i + 1 != len(windows) and i != 0:
                next_window = windows[i + 1]
                self.constrain_approximately_equal(window.A, next_window.A, breakpoint_penalty)
                self.constrain_approximately_equal(window.B, next_window.B, breakpoint_penalty)

            else:
                self.constrain_approximately_equal(window.A, default_ploidy, end_penalty)
                self.constrain_approximately_equal(window.B, default_ploidy, end_penalty)

def reject_outliers(data, m = 2.):
    """http://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list"""
    if len(data) > 1:
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d / mdev if mdev else 0.
        return data[s < m]
    else:
        return data


def get_id():
    """
    Return a unique integer ID (for this execution). Use this to ensure your LP
    variable names are unique.
    
    """
    
    if not hasattr(get_id, "last_id"):
        # Start our IDs at 0, which means the previous ID was -1. Static
        # variables in Python are hard.
        setattr(get_id, "last_id", -1)
        
    # Advance the ID and return the fresh one.
    get_id.last_id += 1
    return get_id.last_id

def plot_it(x, A, B, pngpath, samplename):
    """Plot the inferred copy numbers to pngpath
    x = list of positions"""
    x = np.array(x, dtype="int")
    A = np.array(A, dtype="int")
    B = np.array(B, dtype="int")
    fig = plt.figure()
    plt.axis([x[0], x[-1], 0, 8])
    ax = plt.gca()
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks(range(1,7))
    plt.figtext(0.12, 0.07, "Exon 2")
    plt.figtext(0.81, 0.04, "Intron 2\n(18kb)")
    plt.fill_between(x,A+B, color="blue", alpha=0.7)
    plt.fill_between(x,B, color="red", alpha=0.7)
    plt.scatter(x, A+B, color="blue", label="NOTCH2NL-A")
    plt.scatter(x, B, color="red", label="NOTCH2NL-B")
    r = mpatches.Patch(color="red", label="NOTCH2NL-B", alpha=0.7)
    b = mpatches.Patch(color="blue", label="NOTCH2NL-A", alpha=0.7)
    plt.legend([b, r],["NOTCH2NL-A","NOTCH2NL-B"])
    plt.suptitle("{} NOTCH2NL A/B Copy Number Near Exon 2".format(samplename))
    plt.ylabel("Inferred Copy Number")
    plt.savefig(pngpath, format="png")

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--A", type=str, help="A VCF file", required=True)
    parser.add_argument("--B", type=str, help="B VCF file", required=True)
    parser.add_argument("--N", type=str, help="N2 VCF file", required=True)
    parser.add_argument("--N_positions", type=argparse.FileType("r"), default="n2_normalizing_positions.txt")
    parser.add_argument("--whitelist", type=str, help="whitelist file with probe weights", default="whitelist.txt")
    parser.add_argument("--window", type=int, help="window size to use. Default = 7500bp", default=7500)
    parser.add_argument("--step", type=int, help="step size to use. Default = 1200bp", default=1200)
    parser.add_argument("--png", type=str, help="PNG to write out to", required=True)
    parser.add_argument("--name", type=str, help="Patient Name", required=True)
    parser.add_argument("--breakpoint_penalty", type=float, help="Breakpoint penalty (how easily do we want to open a gene conversion breakpoint?)", default=4)
    parser.add_argument("--data_penalty", type=float, help="data penalty (how accurate are the data?)", default=8)
    parser.add_argument("--end_penalty", type=float, help="end penalty (how often do we expect a full deletion?)", default=3)
    parser.add_argument("--normalize", help="normalize the probe intensities?", action="store_true")
    parser.add_argument("--save_lp_results", help="Save LP results (debugging purposes)?", action="store_true")
    return parser.parse_args()

def main(args):
    args = parse_args(args)

    #map whitelist positions to paralog and weight (and CHM1 pos)
    wl = [x.split() for x in open(args.whitelist) if not x.startswith("#")]
    wl = {x[0]:(x[1],x[2],x[3]) for x in wl}

    #find the Notch2 Adjustment Value based on previously selected positions not in wl
    n2_normalizing = map(int, [x.rstrip() for x in args.N_positions])
    v = [float(x.INFO["ALTFRAC"][0]) for x in vcf.Reader(file(args.N)) if x.POS in n2_normalizing]
    n_adjust = np.mean(reject_outliers(np.array(v))) / 0.8

    
    #make dict mapping positions to adjusted alelle fractions and positions
    value_dict = dict()
    for v in [args.A, args.B]:
        for record in vcf.Reader(file(v)):
            if str(record.POS) in wl:
                paralog, weight, chm1_pos = wl[str(record.POS)]
                if args.normalize is True:
                    val = float(weight) * float(record.INFO["ALTFRAC"][0])
                    value_dict[record.POS] = paralog, val * n_adjust
                else:
                    val = float(record.INFO["ALTFRAC"][0])
                    value_dict[record.POS] = paralog, val * n_adjust

    #initialize problem and model
    problem = SequenceGraphLpProblem()
    model = Model(problem)

    #start populating the model with Windows
    Avals, Bvals = list(), list()
    for x in xrange(region[0], region[1] - args.window, args.step):
        start, end = x, x + args.window
        #make list of values for each paralog
        A = list(); B = list()
        #make list of positions within the current window in the data
        positions = [x for x in value_dict.keys() if x >= start and x < end]
        #iterate over these and assign the values to A or B
        for position in positions:
            paralog, value = value_dict[position]
            if paralog == "A":
                #we multiply by 10 to put it on the same scale as the integer solution
                #I.E. if the value is ~.2 then we are looking for a copy number of 2
                A.append(10 * value)
            elif paralog == "B":
                B.append(10 * value)
        Avals.append(A); Bvals.append(B)

    #now we invert the A and B values because they are on the negative strand
    Avals = Avals[::-1]; Bvals = Bvals[::-1]
    for A, B, pos in zip(Avals, Bvals, range(len(Avals))):
        if len(A) > 0 and len(B) > 0:
            model.add_window(Window(A, B, pos))

    #build and solve the model
    model.sort_windows()
    model.build_model(breakpoint_penalty=args.breakpoint_penalty, data_penalty=args.data_penalty, end_penalty=args.end_penalty)
    problem.solve()

    #find the values for the copy number in each window
    x, A, B = list(), list(), list()
    for pos, window in model.get_windows():
        x.append(pos)
        a, b = window.get_copy_number()
        A.append(a); B.append(b)

    if args.save_lp_results is True:
        #debugging
        tmp = open(os.path.join(os.path.dirname(args.png), args.name + "_ILP_debugging.txt"), "w")
        tmp.write("pos\t(A,B)_data\t(A,B)_inferred\n")
        for pos, window in model.get_windows():
            a,b = window.get_copy_number()
            old_a, old_b = window.get_best_estimate()
            old_a = round(old_a, 3); old_b = round(old_b, 3)
            tmp.write("{}\t{},{}\t{},{}\n".format(pos, old_a, old_b, a, b))
        tmp.close()

    #graph this
    plot_it(x, A, B, args.png, args.name)


if __name__ == "__main__" :
    sys.exit(main(sys.argv))