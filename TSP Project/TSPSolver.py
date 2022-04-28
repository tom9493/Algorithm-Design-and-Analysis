#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq as hq
import copy
import itertools


class State:  # Each state takes up O(n^2) space for the adj matrix
    def __init__(self, matrix, index, currentPath, depth=0, bound=0):
        self.matrix = matrix
        self.path = currentPath
        self.cityIndex = index
        self.depth = depth
        self.bound = bound
        self.key = [-depth, bound, index]

    def reduceCost(self, n):  # Reduce cost algorithm O(n^2)
        low = np.inf
        zero = False
        # Make sure all rows have a zero value, O(n^2)
        for i in range(n):
            zero = False
            low = np.inf
            for j in range(n):
                if self.matrix[i][j] == 0: zero = True
                if self.matrix[i][j] < low: low = self.matrix[i][j]
            if zero:
                continue
            elif low != math.inf:
                for j in range(n):
                    self.matrix[i][j] -= low
                self.bound += low

        # Make sure all columns have a zero value, O(n^2)
        for j in range(n):
            zero = False
            low = np.inf
            for i in range(n):
                if self.matrix[i][j] == 0: zero = True
                if self.matrix[i][j] < low: low = self.matrix[i][j]
            if zero:
                continue
            elif low != math.inf:
                for i in range(n):
                    self.matrix[i][j] -= low
                self.bound += low
        self.updateKey()

    def updateKey(self):  # O(1) time and space
        self.key = [-self.depth, self.bound, self.cityIndex]


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution, 
        time spent to find solution, number of permutations tried during search, the 
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None

        return results

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints: 
        max queue size, total number of states created, and number of pruned states.</returns> 
    '''

    def getInitialMatrix(self, cities, n):  # O(n^2) time and space
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                matrix[i][j] = cities[i].costTo(cities[j])
        return matrix

    def expandState(self, s):
        # If the state only needs to go back to the start city
        if len(s.path) == self.ncities:
            s.bound += s.matrix[s.cityIndex][0]
            if s.bound != math.inf:
                self.stateCount += 1
                self.numSolutions += 1
            if s.bound < self.bssf:
                self.bssf = s.bound
                self.updateQueue()
                self.bestState = s
            return

        # Otherwise expand all other possible routes
        for i in range(1, self.ncities):
            if s.matrix[s.cityIndex][i] != np.inf:
                # Make new state with new bound, O(1) time O(n^2) space
                ns = copy.deepcopy(s)
                ns.depth = s.depth + 1
                ns.cityIndex = i
                ns.bound += ns.matrix[s.cityIndex][i]
                ns.path.append(self.cities[i])
                self.stateCount += 1

                # Infinite out appropriate columns
                for j in range(self.ncities):
                    ns.matrix[s.cityIndex][j] = math.inf
                for j in range(self.ncities):
                    ns.matrix[j][i] = math.inf
                ns.matrix[i][s.cityIndex] = math.inf

                ns.reduceCost(self.ncities)
                if ns.bound <= self.bssf:
                    hq.heappush(self.Q, (ns.key, ns))
                else:
                    self.pruneCount += 1

        if len(self.Q) > self.maxQueueSize: self.maxQueueSize = len(self.Q)

    def updateQueue(self):
        r = 0
        for i in range(len(self.Q)):
            if self.Q[i][1].bound >= self.bssf:
                r += 1
        for i in range(len(self.Q) - r):
            if self.Q[i][1].bound >= self.bssf:
                self.Q[i] = self.Q[-1]
                hq.heappop(self.Q)
                if i < len(self.Q):
                    hq._siftup(self.Q, i)
                    hq._siftdown(self.Q, 0, i)
                self.pruneCount += 1

    def branchAndBound(self, time_allowance=60.0):
        self.cities = self._scenario.getCities()
        self.ncities = len(self.cities)
        path = []
        path.append(self.cities[0])
        initialState = State(self.getInitialMatrix(self.cities, self.ncities), 0, path)
        initialState.reduceCost(self.ncities)
        default = self.defaultRandomTour()

        self.bestState = None
        self.numSolutions = 0
        self.stateCount = 0
        self.pruneCount = 0
        self.maxQueueSize = 1
        self.bssf = default['cost']
        self.Q = []

        hq.heappush(self.Q, (initialState.key, initialState))
        start_time = time.time()
        # Branch and Bound algorithm
        times = 0
        while len(self.Q) != 0 and time.time() - start_time < time_allowance:
            times += 1
            current = hq.heappop(self.Q)[1]
            self.expandState(current)
        end_time = time.time()

        print(times)
        results = {}
        results['cost'] = self.bssf
        results['time'] = end_time - start_time
        results['count'] = self.numSolutions
        results['max'] = self.maxQueueSize
        results['total'] = self.stateCount
        results['pruned'] = self.pruneCount
        if self.bestState != None: results['soln'] = TSPSolution(self.bestState.path)

        else:
            results['soln'] = default['soln']
        return results

    ''' <summary>
            This is the entry point for the greedy solver, which you must implement for 
            the group project (but it is probably a good idea to just do it for the branch-and
            bound project as a way to get your feet wet).  Note this could be used to find your
            initial BSSF.
            </summary>
            <returns>results dictionary for GUI that contains three ints: cost of best solution, 
            time spent to find best solution, total number of solutions found, the best
            solution found, and three null values for fields not used for this 
            algorithm</returns> 
        '''

    def greedy(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found during search, the 
        best solution found.  You may use the other three field however you like.
        algorithm</returns> 
    '''

    def getProbability(self, i, j):     # O(n) time, O(1) space
        tij = self.phMtx[i][j]  # Intensity of pheromone trail from i to j
        if self.adjMtx[i][j] != 0: nij = 1 / self.adjMtx[i][j]  # Visibility from i to j ( 1/distance )
        else: return math.inf
        sumAllowed = 0
        for k in range(self.ncities):   # O(n) iterations, constant work
            if self.adjMtx[i][k] != math.inf:
                tik = self.phMtx[i][k]
                if self.adjMtx[i][k] != 0: nik = 1 / self.adjMtx[i][k]
                else: nik = 1 / 0.00001
                sumAllowed += tik ** self.alpha * nik ** self.beta
        return (tij ** self.alpha * nij ** self.beta) / sumAllowed

    def updatePheromone(self):  # O(n^2) time, O(1) space
        self.phChgMtxTotal = self.phChgMtxTotal + self.phChgMtx
        self.phChgMtx = np.zeros((self.ncities, self.ncities))
        for i in range(self.ncities):
            for j in range(self.ncities):
                self.phMtx[i][j] = self.p * self.phMtx[i][j] + self.phChgMtxTotal[i][j]

    def infOut(self, i, j, mtx):       # O(n) time, O(1) space
        for k in range(self.ncities):
            mtx[i][k] = math.inf
        for k in range(self.ncities):
            mtx[k][j] = math.inf
        mtx[j][i] = math.inf
        if i == self.startCityIndex:
            for k in range(self.ncities):
                mtx[k][i] = math.inf
        return mtx

    def makePhMtx(self, path, totalDist):   # O(n) time, O(1) space
        for i in range(self.ncities):
            if i == self.ncities - 1:
                a = path[i][0]
                b = path[0][0]
            else:
                a = path[i][0]
                b = path[i+1][0]
            self.phChgMtx[a][b] = 1 / totalDist

    def fancy(self, time_allowance=60.0):   # This implementation is O(n^4) time, O(n^2) space, moderately optimal
        self.cities = self._scenario.getCities()     # All these O(n^2) space
        self.ncities = len(self.cities)
        self.adjMtx = self.getInitialMatrix(self.cities, self.ncities)
        self.phMtx = np.full((self.ncities, self.ncities), 1.0)
        self.phChgMtxTotal = np.zeros((self.ncities, self.ncities))
        self.phChgMtx = np.zeros((self.ncities, self.ncities))
        self.alpha = 1  # Influence of pheromone on probability
        self.beta = 5   # Influence of distance on probability
        self.p = 0.5    # Evaporation factor of pheromone

        iterations = 0
        bestPathLength = math.inf
        bestPath = None
        start_time = time.time()             # constant factor iterations
        while iterations < 10 and time.time() - start_time < time_allowance: # O(n) iterations
            self.startCityIndex = random.randrange(self.ncities)
            path = []
            path.append([self.startCityIndex, self.cities[self.startCityIndex]])
            cityIndex = self.startCityIndex
            totalDist = 0
            mtx = self.adjMtx.copy()
            for i in range(self.ncities-1):     # O(n) iterations, so O(n^3)
                bestJ = None
                bestP = 0
                for j in range(self.ncities):   # O(n) iterations
                    if mtx[cityIndex][j] != math.inf:
                        probability = self.getProbability(cityIndex, j) # O(n) time
                        if probability > bestP:
                            bestP = probability
                            bestJ = j
                if cityIndex != bestJ and bestJ != None:
                    totalDist = totalDist + mtx[cityIndex][bestJ]
                    if i == self.ncities - 2: totalDist = totalDist + self.adjMtx[bestJ][self.startCityIndex]
                    mtx = self.infOut(cityIndex, bestJ, mtx) # O(n) time
                    cityIndex = bestJ
                    path.append([cityIndex, self.cities[cityIndex]])
                else: totalDist = math.inf

            if totalDist != math.inf:
                self.makePhMtx(path, totalDist)     # O(n) time
                self.updatePheromone()              # O(n^2) time
            if totalDist < bestPathLength:
                bestPathLength = totalDist
                bestPath = path
            elif totalDist == math.inf: continue
            iterations += 1

        end_time = time.time()
        solution = []
        for i in range(self.ncities):       # O(n) time and space
            solution.append(bestPath[i][1])

        results = {}
        results['soln'] = TSPSolution(solution)
        results['cost'] = bestPathLength
        results['time'] = end_time - start_time
        results['count'] = iterations
        results['max'] = None
        results['total'] = None
        results['pruned'] = None

        return results


