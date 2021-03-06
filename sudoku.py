## Solve Every Sudoku Puzzle

## See http://norvig.com/sudoku.html

## Throughout this program we have:
##   r is a row,    e.g. 'A'
##   c is a column, e.g. '3'
##   s is a square, e.g. 'A3'
##   d is a digit,  e.g. '9'
##   u is a unit,   e.g. ['A1','B1','C1','D1','E1','F1','G1','H1','I1']
##   grid is a grid,e.g. 81 non-blank chars, e.g. starting with '.18...7...
##   values is a dict of possible values, e.g. {'A1':'12349', 'A2':'8', ...}

import sys
import copy
import random as rand
from math import exp

attempt_cnt = 0  # computeur du nombre de tentatives

def cross(A, B):
    "Cross product of elements in A and elements in B."
    return [a+b for a in A for b in B]

digits   = '123456789'
rows     = 'ABCDEFGHI'
cols     = digits
squares  = cross(rows, cols)
unitlist = ([cross(rows, c) for c in cols] +
            [cross(r, cols) for r in rows] +
            [cross(rs, cs) for rs in ('ABC','DEF','GHI') for cs in ('123','456','789')])
units = dict((s, [u for u in unitlist if s in u])
             for s in squares)
peers = dict((s, set(sum(units[s],[]))-set([s]))
             for s in squares)

col_units = unitlist[:9]
row_units = unitlist[9:18]
sqr_units = unitlist[18:]

################ Unit Tests ################

def test():
    "A set of tests that must pass."
    assert len(squares) == 81
    assert len(unitlist) == 27
    assert all(len(units[s]) == 3 for s in squares)
    assert all(len(peers[s]) == 20 for s in squares)
    assert units['C2'] == [['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2'],
                           ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'],
                           ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']]
    assert peers['C2'] == set(['A2', 'B2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2',
                               'C1', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                               'A1', 'A3', 'B1', 'B3'])
    print 'All tests pass.'



############################ Hill Climbing ###################################################

def hill_climbing(values):
    global attempt_cnt
    attempt_cnt = 0

    progress = True

    square_units = [cross(rs, cs) for rs in ('ABC','DEF','GHI') for cs in ('123','456','789')]
    empty_squares = [s for s in squares if values[s] in '0.']

    ## For each empty square in a square unit, assign a digit not already in this unit.
    for u in square_units:
        ds = set(digits) - set(values[s] for s in u)
        for s in u:
            if values[s] in '0.':
                values[s] = random.sample(ds, 1)[0]
                ds.remove(values[s])

    while progress:
        conflicts = get_conflicts(values)
        initial_best = len(conflicts)
        best = initial_best

        rand.shuffle(square_units)
        progress = False

        for u in square_units:
            prospect = set()

            for s1 in u:
                for s2 in u:
                    if s1 != s2 and s1 in empty_squares and s2 in empty_squares:
                        new_values = copy.deepcopy(values)
                        new_values[s1], new_values[s2] = new_values[s2], new_values[s1]
                        attempt_cnt += 1

                        new_conflicts = get_conflicts(new_values)
                        if len(new_conflicts) < best:
                            prospect = set()
                            prospect.add((s1, s2))
                            best = len(new_conflicts)
                            # print "best: " + str(best)
                        elif len(new_conflicts) == best:
                            prospect.add((s1, s2))

            # if there are some prospects, swap the squares in a random prospect
            if len(prospect) > 0:
                s1, s2 = prospect.pop()
                values[s1], values[s2] = values[s2], values[s1]

            if best == 0:
                print "win!!!!!!!!!!!!!!!"
                return values  ## Solved!
            elif best != initial_best:
                progress = True

    print "----------- fail: " + str(best)
    return False

############################### Simulated Annealing ###################################################

# Simulated annealing utilisant hill_climbing
def simulated_annealing(values):
    global attempt_cnt
    attempt_cnt = 0

    best = 0

    progress = True

    square_units = [cross(rs, cs) for rs in ('ABC', 'DEF', 'GHI') for cs in ('123', '456', '789')]
    empty_squares = [s for s in squares if values[s] in '0.']

    ## For each empty square in a square unit, assign a digit not already in this unit.
    for u in square_units:
        ds = set(digits) - set(values[s] for s in u)
        for s in u:
            if values[s] in '0.':
                values[s] = random.sample(ds, 1)[0]
                ds.remove(values[s])

    T = 2.0

    while progress:
        conflicts = get_conflicts(values)
        initial_best = len(conflicts)
        best = initial_best

        rand.shuffle(square_units)
        progress = False

        for u in square_units:
            prospect = set()

            for s1 in u:
                for s2 in u:
                    if s1 != s2 and s1 in empty_squares and s2 in empty_squares:
                        new_values = copy.deepcopy(values)
                        new_values[s1], new_values[s2] = new_values[s2], new_values[s1]
                        attempt_cnt += 1

                        new_conflicts = get_conflicts(new_values)
                        if len(new_conflicts) < best:
                            prospect = set()
                            prospect.add((s1, s2))
                            best = len(new_conflicts)
                            # print "best: " + str(best)
                        elif len(new_conflicts) == best:
                            prospect.add((s1, s2))

            # if there are some prospects, swap the squares in a random prospect
            if len(prospect) > 0:
                s1, s2 = prospect.pop()
                values[s1], values[s2] = values[s2], values[s1]

            if best == 0:
                print "win!!!!!!!!!!!!!!!"
                return values  ## Solved!
            elif best != initial_best:
                progress = True

        # find a square unit in which there's enough non-fixed squares for a swap
        swap_prospects = set()
        while len(swap_prospects) < 2:
            u = random.choice(sqr_units)
            swap_prospects = set.intersection(set(u), empty_squares)

        to_swap = random.sample(swap_prospects, 2)
        s1, s2 = to_swap[0], to_swap[1]

        # swap the selected squares
        new_values = copy.deepcopy(values)
        new_values[s1], new_values[s2] = new_values[s2], new_values[s1]

        print "best: " + str(best) + "   T: " + str(T)

        if not progress:
            new_conflicts = get_conflicts(new_values)
            dE = best - len(new_conflicts)

            if True or dE > 0 or exp(dE / T) > rand.uniform(0, 1):
                progress = True
                values[s1], values[s2] = values[s2], values[s1]
                best -= dE

        T *= 0.99

    print "----------- fail: " + str(best)
    return False


def get_conflicts(values):
    """Return a list of squares causing conflict."""
    return [s for s in squares if values[s] in [values[s2] for s2 in peers[s]]]


############### Heuristique 1 #####################

# Simulated annealing inspire de Lewis, heuristique basee sur le nombre de chiffres manquants dans les unites
def sa_heur1(initial_values):
    global attempt_cnt
    attempt_cnt = 0

    empty_squares = set([s for s in squares if initial_values[s] in '0.'])

    # will loop until the solution is found
    while True:
        # For each empty square in a square unit, assign a digit not already in this unit.
        values = copy.deepcopy(initial_values)
        for u in sqr_units:
            ds = set(digits) - set(values[s] for s in u)
            for s in u:
                if values[s] in '0.':
                    values[s] = random.sample(ds, 1)[0]
                    ds.remove(values[s])

        # get the count of missing numbers in each col/row and sum them to get the cost
        cols_cost = [9 - len(set(values[s] for s in u)) for u in col_units]
        rows_cost = [9 - len(set(values[s] for s in u)) for u in row_units]
        current_cost = sum(cols_cost) + sum(rows_cost)

        T = 3.0  # temperature

        while T > 0.1:
            # find a square unit in which there's enough non-fixed squares for a swap
            swap_prospects = set()
            while len(swap_prospects) < 2:
                u = random.choice(sqr_units)
                swap_prospects = set.intersection(set(u), empty_squares)

            # swap two squares in selected unit
            s1, s2 = random.sample(swap_prospects, 2)
            values[s1], values[s2] = values[s2], values[s1]
            attempt_cnt += 1

            # keep a copy of cost arrays in case we do not swap
            rows_cost_cpy = copy.deepcopy(rows_cost)
            cols_cost_cpy = copy.deepcopy(cols_cost)

            # get the cost of the new position
            new_cost = compute_cost_h1(values, rows_cost, cols_cost, s1, s2)

            # check if we accept the swap or reverse it
            dE = current_cost - new_cost
            if dE > 0 or exp(dE / T) > rand.uniform(0, 1):
                current_cost = new_cost
            else:
                values[s1], values[s2] = values[s2], values[s1]
                rows_cost = rows_cost_cpy
                cols_cost = cols_cost_cpy

            # update temperature if solution not found, otherwise return solution
            if current_cost != 0:
                T *= 0.9999
            else:
                print "win!!!!!!!!!!!!!!!"
                return values

            # print "cost: " + str(current_cost)

        print "reheat   ---   attempts: " + str(attempt_cnt) + "   cost: " + str(current_cost)


def compute_cost_h1(values, rows_cost, cols_cost, s1, s2):
    affected_rows = set([ord(s1[0]), ord(s2[0])])
    affected_cols = set([ord(s1[1]), ord(s2[1])])

    for r in affected_rows:
        rows_cost[r - 65] = 9 - len(set([values[s] for s in row_units[r - 65]]))

    for c in affected_cols:
        cols_cost[c - 49] = 9 - len(set([values[s] for s in col_units[c - 49]]))

    return sum(rows_cost) + sum(cols_cost)


############### Heuristique 2 ####################

# Simulated annealing avec heuristique basee sur le nombre de cases en conflit
def sa_heur2(values):
    global attempt_cnt
    attempt_cnt = 0

    square_units = [cross(rs, cs) for rs in ('ABC','DEF','GHI') for cs in ('123','456','789')]
    empty_squares = [s for s in squares if values[s] in '0.']

    T = 0.1  # temperature

    ## For each empty square in a square unit, assign a digit not already in this unit.
    for u in square_units:
        ds = set(digits) - set(values[s] for s in u)
        for s in u:
            if values[s] in '0.':
                values[s] = random.sample(ds, 1)[0]
                ds.remove(values[s])

    while T > 0.0001:
        conflicts = get_conflicts(values)
        best = len(conflicts)

        # find a square unit in which there's enough non-fixed squares for a swap
        swap_prospects = set()
        while len(swap_prospects) < 2:
            u = random.choice(sqr_units)
            swap_prospects = set.intersection(set(u), empty_squares)

        # choose two squares to swap in selected square unit
        to_swap = random.sample(swap_prospects, 2)
        s1, s2 = to_swap[0], to_swap[1]

        # swap the selected squares
        new_values = copy.deepcopy(values)
        new_values[s1], new_values[s2] = new_values[s2], new_values[s1]
        attempt_cnt += 1

        new_conflicts = get_conflicts(new_values)
        dE = best - len(new_conflicts)

        if dE > 0 or exp(dE/T) > rand.uniform(0, 1):
            values[s1], values[s2] = values[s2], values[s1]
            best -= dE

        # print "best: " + str(best)
        # print "T: " + str(T)

        if best == 0:
            print "win!!!!!!!!!!!!!!!"
            return values  ## Solved!
        else:
            T *= 0.999

    print "----------- fail: " + str(best)
    return False


################ Parse a Grid ################

def parse_grid(grid):
    """Convert grid to a dict of possible values, {square: digits}, or
    return False if a contradiction is detected."""
    ## To start, every square can be any digit; then assign values from the grid.
    values = dict((s, digits) for s in squares)
    for s,d in grid_values(grid).items():
        if d in digits and not assign(values, s, d):
            return False ## (Fail if we can't assign d to square s.)
    return values

def grid_values(grid):
    "Convert grid into a dict of {square: char} with '0' or '.' for empties."
    chars = [c for c in grid if c in digits or c in '0.']
    assert len(chars) == 81
    return dict(zip(squares, chars))

################ Constraint Propagation ################

def assign(values, s, d):
    """Eliminate all the other values (except d) from values[s] and propagate.
    Return values, except return False if a contradiction is detected."""
    other_values = values[s].replace(d, '')
    if all(eliminate(values, s, d2) for d2 in other_values):
        global attempt_cnt
        attempt_cnt += 1
        return values
    else:
        return False

def eliminate(values, s, d):
    """Eliminate d from values[s]; propagate when values or places <= 2.
    Return values, except return False if a contradiction is detected."""
    if d not in values[s]:
        return values ## Already eliminated
    values[s] = values[s].replace(d,'')
    ## (1) If a square s is reduced to one value d2, then eliminate d2 from the peers.
    if len(values[s]) == 0:
        return False ## Contradiction: removed last value
    elif len(values[s]) == 1:
        d2 = values[s]
        if not all(eliminate(values, s2, d2) for s2 in peers[s]):
            return False
    ## (2) If a unit u is reduced to only one place for a value d, then put it there.
    for u in units[s]:
        dplaces = [s for s in u if d in values[s]]
        if len(dplaces) == 0:
            return False ## Contradiction: no place for this value
        elif len(dplaces) == 1:
            # d can only be in one place in unit; assign it there
            if not assign(values, dplaces[0], d):
                return False
    return values

################ Display as 2-D grid ################


def display(values):
    "Display these values as a 2-D grid."
    width = 1+max(len(values[s]) for s in squares)
    line = '+'.join(['-'*(width*3)]*3)
    for r in rows:
        print ''.join(values[r+c].center(width)+('|' if c in '36' else '')
                      for c in cols)
        if r in 'CF': print line
    print

################ Search ################

def solve(grid, method): 
    if method == 'norvig':
        return search(parse_grid(grid))
    elif method == 'hc':
        return hill_climbing(grid_values(grid))
    elif method == 'sa':
        return simulated_annealing(grid_values(grid))
    elif method == 'sa_heur1':
        return sa_heur1(grid_values(grid))
    elif method == 'sa_heur2':
        return sa_heur2(grid_values(grid))
    else:
        print "\nWRONG COMMAND LINE ARGUMENT\n"

def search(values):
    "Using depth-first search and propagation, try all possible values."
    if values is False:
        return False ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values ## Solved!
    ## Chose the unfilled square s with the fewest possibilities
    n,s = min((len(values[s]), s) for s in squares if len(values[s]) > 1)
    return some(search(assign(values.copy(), s, d))
                for d in values[s])

################ Utilities ################


def some(seq):
    "Return some element of seq that is true."
    for e in seq:
        if e: return e
    return False


def from_file(filename, sep='\n'):
    "Parse a file into a list of strings, separated by sep."
    return file(filename).read().strip().split(sep)


def shuffled(seq):
    "Return a randomly shuffled copy of the input sequence."
    seq = list(seq)
    random.shuffle(seq)
    return seq

################ System test ################


import time, random


def solve_all(grids, method, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results.
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""
    def time_solve(grid):
        global attempt_cnt
        start = time.clock()
        attempt_cnt = 0
        values = solve(grid, method)
        print attempt_cnt
        t = time.clock()-start
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print '(%.2f seconds)\n' % t
        return (t, solved(values))
    times, results = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N > 1:
        print "Solved %d of %d %s puzzles (avg %.2f secs (%d Hz), max %.2f secs)." % (
            sum(results), N, name, sum(times)/N, N/sum(times), max(times))

def solved(values):
    "A puzzle is solved if each unit is a permutation of the digits 1 to 9."
    def unitsolved(unit): return set(values[s] for s in unit) == set(digits)
    return values is not False and all(unitsolved(unit) for unit in unitlist)

def random_puzzle(N=17):
    """Make a random puzzle with N or more assignments. Restart on contradictions.
    Note the resulting puzzle is not guaranteed to be solvable, but empirically
    about 99.8% of them are solvable. Some have multiple solutions."""
    values = dict((s, digits) for s in squares)
    for s in shuffled(squares):
        if not assign(values, s, random.choice(values[s])):
            break
        ds = [values[s] for s in squares if len(values[s]) == 1]
        if len(ds) >= N and len(set(ds)) >= 8:
            return ''.join(values[s] if len(values[s])==1 else '.' for s in squares)
    return random_puzzle(N) ## Give up and make a new puzzle

grid1  = '003020600900305001001806400008102900700000008006708200002609500800203009005010300'
grid2  = '4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......'
hard1  = '.....6....59.....82....8....45........3........6..3.54...325..6..................'
    
if __name__ == '__main__':
    test()
    assert len(sys.argv) != 1
    # solve_all(from_file("easy50.txt", '========'), sys.argv[1], "easy", None)
    # solve_all(from_file("top95.txt"), sys.argv[1], "hard", None)
    # solve_all(from_file("1000sudoku.txt"), sys.argv[1], "hard", None)
    # solve_all(from_file("hardest.txt"), sys.argv[1], "hardest", None)
    # solve_all([random_puzzle() for _ in range(99)], "random", 100.0)
    solve_all(from_file("100sudoku.txt"), sys.argv[1], "--- 100 sudoku ---", None)



## References used:
## http://www.scanraid.com/BasicStrategies.htm
## http://www.sudokudragon.com/sudokustrategy.htm
## http://www.krazydad.com/blog/2005/09/29/an-index-of-sudoku-strategies/
## http://www2.warwick.ac.uk/fac/sci/moac/currentstudents/peter_cock/python/sudoku/
