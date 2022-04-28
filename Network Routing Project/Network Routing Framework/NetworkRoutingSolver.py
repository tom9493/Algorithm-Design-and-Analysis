#!/usr/bin/python3
import math

from CS312Graph import *
import time


class NetworkRoutingSolver:
    def __init__(self):
        pass

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network

    def getShortestPath(self, destIndex):
        self.dest = destIndex
        path_edges = []
        total_length = 0
        currentNode = self.network.nodes[destIndex]
        prevNode = self.network.nodes[self.prev[currentNode.node_id]]
        while currentNode.node_id != self.source:
            for i in range(3):
                if prevNode.neighbors[i].dest == currentNode:
                    edge = prevNode.neighbors[i]
                    path_edges.append((edge.src.loc, edge.dest.loc, '{:.0f}'.format(edge.length)))
                    total_length += edge.length
                    currentNode = prevNode
                    prevNode = self.network.nodes[self.prev[currentNode.node_id]]
        return {'cost': total_length, 'path': path_edges}

    def computeShortestPaths(self, srcIndex, use_heap=False):
        self.dist = []
        self.prev = []
        self.source = srcIndex
        self.use_heap = use_heap
        t1 = time.time()

        for i in range(len(self.network.nodes)): # O(V)
            self.dist.append(math.inf)
            self.prev.append(-1)
            if i == srcIndex:
                self.dist[i] = 0

        H = self.makeQueue()         # Heap, array: O(V)
        while self.queueSize > 0:
            u = self.network.nodes[self.deleteMin(H)]   # Heap: O(logV), Array: O(V)
            for i in range(3):
                v = u.neighbors[i].dest
                if self.dist[v.node_id] > self.dist[u.node_id] + u.neighbors[i].length:
                    self.dist[v.node_id] = self.dist[u.node_id] + u.neighbors[i].length
                    self.prev[v.node_id] = u.node_id
                    self.decreaseKey(H, v.node_id)      # Heap: O(logV), Array: O(1)

        t2 = time.time()
        return (t2 - t1)

    # Makes the priority queue, implementing the heap or array depending on use_heap
    def makeQueue(self):
        if not self.use_heap:                 # Array: O(V)
            H = []
            self.queueSize = len(self.dist)
            for i in range(len(self.dist)):
                H.append(-1)
                if self.dist[i] == 0:
                    H[i] = 0

        else:                                 # Heap: O(V)
            self.indices = []
            for i in range(len(self.dist)):
                self.indices.append(-1)

            H = []
            H.append([self.source, 0])
            self.indices[self.source] = 0
            self.queueSize = 1

        return H

    # Deletes the smallest distance value in priority queue. Returns its node_id/index
    def deleteMin(self, H):
        if not self.use_heap:  # Array deleteMin O(V)
            min = math.inf
            index = -1
            for i in range(len(H)):
                if H[i] == -1:
                    continue
                if H[i] == 0:
                    H[i] = -1
                    return i
                if H[i] < min:
                    min = H[i]
                    index = i
            H[index] = -1

        else:                    # Heap deleteMin O(logV)
            index = H[0][0]
            H[0] = H[self.queueSize - 1]
            H.pop(self.queueSize - 1)
            if len(H) != 0: self.shiftDown(H)   # shiftDown: O(logV)

        self.queueSize -= 1
        return index

    # If a closer distance has been found, decrease distance key in priority queue
    def decreaseKey(self, H, nodeId):

        if not self.use_heap:   # Array decreaseKey O(1)
            H[nodeId] = self.dist[nodeId]

        else:                   # Heap decreaseKey AND insert O(logV)
            if self.indices[nodeId] == -1:
                self.insert(H, nodeId, self.dist[nodeId])
            else:
                H[self.indices[nodeId]][1] = self.dist[nodeId]
                self.shiftUp(H, self.indices[nodeId])

    # Inserts a node that doesn't exist into the priority queue(heap)
    def insert(self, H, index, dist): # Heap insert O(logV)
        H.append([index, dist])
        self.queueSize += 1
        self.shiftUp(H, self.queueSize - 1) # shiftUp O(logV)

    # These functions allow us to find the indices we need in the heap. All O(1)
    def getParent(self, childIndex):
        return (childIndex - 1) // 2

    def getLeftChild(self, parentIndex):
        return (2 * parentIndex) + 1

    def getRightChild(self, parentIndex):
        return (2 * parentIndex) + 2

    # Function needed for deleteMin. Pop the top node, replace it with the last node
    # in the array we are using as the heap, then shift it down until it is sorted
    def shiftDown(self, H):         # O(logV) worst case
        ni = 0
        li = self.getLeftChild(ni)
        ri = self.getRightChild(ni)
        node = H[ni]

        if li < len(H): l = H[li]   # left and right index assignments
        else: l = [-1, math.inf]    # if the indices don't exist in the
        if ri < len(H): r = H[ri]   # priority queue, set to inf
        else: r = [-1, math.inf]

        while node[1] > l[1] or node[1] > r[1]:
            if l[1] < r[1]:
                H[ni], H[li] = H[li], H[ni]
                self.indices[node[0]] = li
                self.indices[l[0]] = ni
                ni = li
            else:
                H[ni], H[ri] = H[ri], H[ni]
                self.indices[node[0]] = ri
                self.indices[r[0]] = ni
                ni = ri
            li = self.getLeftChild(ni)
            ri = self.getRightChild(ni)
            node = H[ni]
            if li < len(H): l = H[li]
            else: l = [-1, math.inf]
            if ri < len(H): r = H[ri]
            else: r = [-1, math.inf]

    # Function needed for insert or decreaseKey. Takes given index of a node that
    # has been adjusted and shifts it up according to its smaller distance value
    def shiftUp(self, H, index):      # O(logV) worst case
        ni = index
        pi = self.getParent(ni)
        node = H[ni]
        p = H[pi]
        while node[1] < p[1]:
            H[ni], H[pi] = H[pi], H[ni]
            self.indices[node[0]] = pi
            self.indices[p[0]] = ni
            ni = pi
            pi = self.getParent(ni)
            node = H[ni]
            p = H[pi]

