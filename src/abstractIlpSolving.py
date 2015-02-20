import pulp, os

"""
All code in this file written by Adam Novak
"""


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
        self.is_solved = None

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

    def is_solved(self):
        return self.is_solved

    def solve(self, save=None):
        """
        Solve the LP problem with CBC
        
        If save is specified, it is a filename to which to save the LP problem
        in LP format.
        
        You may only solve a SequenceGraphLpProblem once.
        
        """
        # Set up the penalties described by the penalty tree
        self.penalties.set_objective(self.problem)
        # Solve the problem
        status = self.problem.solve(pulp.COIN_CMD(path=os.path.join(os.path.expanduser("~"), "bin", "cbc")))
        if status != pulp.constants.LpStatusOptimal:
            self.is_solved = False
            raise Exception("Unable to solve problem optimally.")
        else:
            self.is_solved = True
