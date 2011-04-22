'''
Created on 01-01-2011

@author: lisu
'''
from math import sqrt
import sys
from pydoc import deque
from timeit import itertools

def delLine():
    for _ in range(255):
        sys.stdout.write('\010') #backspace
        
def clearFile(filename):
    f = open(filename, 'w')
    f.close()

def posmax(seq, key=None): 
    """Returns position of the max value on the list"""
    if key is None:
        m = seq[0]
        index = 0
        for i, x in enumerate(seq):
            if x > m:
                m = x
                index = i
        return index
    else:
        return NotImplemented

def posmin(seq, key=None): 
    if key is None:
        m = seq[0]
        index = 0
        for i, x in enumerate(seq):
            if x < m:
                m = x
                index = i
        return index
    else:
        return NotImplemented
    
def concatenate(Lists):
    """ [[1],[2],[3]] --> [1,2,3] """
    resList = []
    for l in Lists:
        resList.extend(l)
    return resList

def bruteForce(listLength, elements):
    """ 
    for 1, [0,1] will return
    [0,0]
    [0,1]
    [1,0]
    [1,1]
    """
    if listLength > 1:
        for element in elements:
            for i in bruteForce(listLength - 1, elements):
                yield [element] + i
    elif listLength == 1:
        for i in elements:
            yield [i]

def distance(x1, x2):
    """Calculates distance between vectors x1 and x2"""
    return sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2)
        #return sqrt(sum([(x1[i] - x2[i])**2 for i in range(len(x1))]))

def vectorLength(x1):
    return distance([0,0,0], x1)

def movingAverage(iterable, n=3):
    # http://en.wikipedia.org/wiki/Moving_average
    it = list(iterable)

    for elem in range(len(it)):
        s = 0
        num = positiveNumOrZero(elem - n)
        windowSize = min(len(it[num:elem]), len(it[elem + 1:]))
        
        for i in range(1, windowSize + 1): 
            s += it[elem + i ] + it[elem - i]
        s += it[elem] 
        yield s / float(2 * windowSize + 1)

def positiveNumOrZero(number):
    if number < 0: 
        return 0
    else: 
        return number   