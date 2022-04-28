from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25


#
# This is the class you have to complete.
#

class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This function implements a divide and conquer algorithm for finding a convex hull of an array of sorted points
    # Assume all Math functions are all constant time (division, subtraction, etc)
    def convexHullSolver(self, a):				# Whole thing O(nlogn)
        if len(a) < 4:                          # Space complexity also O(nlogn)
            return self.makeBaseHull(a)         # makeBaseHull function O(1)
        else:
            halfIndex = len(a) // 2             # Division is O(1)
            lowerArray = a[0:halfIndex]         # Each recursion takes O(n) more space
            upperArray = a[halfIndex: len(a)]
            lowerArray = self.convexHullSolver(lowerArray)          # O(n)
            upperArray = self.convexHullSolver(upperArray)          # O(n)
            newArray = self.mergeHulls(lowerArray, upperArray)      # mergeHull O(n)
        return newArray

    # This takes an array of 3 or fewer points, returning it if it is 1 or two points in size, and returning the
    # clockwise order if it is 3 points in size. This is all in constant time.
    def makeBaseHull(self, a):      # Whole thing O(1)
        if len(a) < 3:
            return a
        else:
            if self.findSlope(a[0], a[1]) > self.findSlope(a[0], a[2]):
                return a
            else:
                return [a[0], a[2], a[1]]

    # Constant time function to find the slope of a line between two given points
    def findSlope(self, point1, point2):        # Whole thing O(1)
        x1, y1 = point1.x(), point1.y()
        x2, y2 = point2.x(), point2.y()
        return (y2-y1)/(x2-x1)

    # This function merges two convex hulls together regardless of its size. It takes the left-most point on the
    # right hull and the right-most point on the left hull and compares slopes incrementally to find the upper and
    # lower tangents. A new array of clockwise points in formed representing the new convex hull. The other points
    # not included in the hull are forgotten, giving us a slight improvement in time.
    def mergeHulls(self, lt, rt):               # Whole thing O(n)
        lIndex = self.getRightPoint(lt)         # getRightPoint function is O(n)
        rIndex = 0
        slope = self.findSlope(lt[lIndex], rt[rIndex])

        # All of these while loops will be less than O(n) with constant work each loop so O(n) and O(1) space complexity
        while (slope > self.findSlope(lt[lIndex - 1], rt[rIndex]) or slope <
               self.findSlope(lt[lIndex], rt[(rIndex + 1) % len(rt)])):
            while slope > self.findSlope(lt[lIndex - 1], rt[rIndex]):
                slope = self.findSlope(lt[lIndex - 1], rt[rIndex])
                lIndex = lIndex - 1
            while slope < self.findSlope(lt[lIndex], rt[(rIndex + 1) % len(rt)]):
                slope = self.findSlope(lt[lIndex], rt[(rIndex + 1) % len(rt)])
                rIndex = (rIndex + 1) % len(rt)

        luTangent = lIndex % len(lt)            # O(1) division
        ruTangent = rIndex % len(rt)            # O(1) division
        lIndex = self.getRightPoint(lt)         # O(n) getRightPoint
        rIndex = 0
        slope = self.findSlope(lt[lIndex], rt[rIndex])

        # All of these while loops will be less than O(n) with constant work each loop so O(n) and O(1) space complexity
        while (slope > self.findSlope(lt[lIndex], rt[rIndex - 1]) or slope <
               self.findSlope(lt[(lIndex + 1) % len(lt)], rt[rIndex])):
            while slope > self.findSlope(lt[lIndex], rt[rIndex - 1]):
                slope = self.findSlope(lt[lIndex], rt[rIndex - 1])
                rIndex = rIndex - 1
            while slope < self.findSlope(lt[(lIndex + 1) % len(lt)], rt[rIndex]):
                slope = self.findSlope(lt[(lIndex + 1) % len(lt)], rt[rIndex])
                lIndex = (lIndex + 1) % len(lt)

        rIndex = rIndex % len(rt)
        newArray = []

        # No more than O(n) for these loops
        for i in range(luTangent + 1):
            newArray.append(lt[i])
        while ruTangent != rIndex:
            newArray.append(rt[ruTangent])
            ruTangent = (ruTangent + 1) % len(rt)
        newArray.append(rt[rIndex])
        while lIndex != 0:
            newArray.append(lt[lIndex])
            lIndex = (lIndex + 1) % len(lt)

        return newArray

    # This finds the point with the right most x-coordinate in an array and returns its index
    def getRightPoint(self, a):             # O(n) time
        index = 0
        for i in range(len(a)):             # Loop makes this O(n)
            if a[i].x() > a[index].x():
                index = i
        return index

    # The given function for the project. Whole thing O(nlogn) time AND space complexity after sort and solve
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        points.sort(key=lambda point: point.x())            # O(nlogn) sort
        t2 = time.time()
        t3 = time.time()

        points = self.convexHullSolver(points)              # O(nlogn) time AND space complexity
        polygon = [QLineF(points[i], points[(i + 1) % len(points)]) for i in range(len(points))]
        t4 = time.time()

        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))
