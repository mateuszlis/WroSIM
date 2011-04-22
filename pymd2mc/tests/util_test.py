'''
Created on 22-01-2011

@author: lisu
'''

from random import randint
import unittest


from utils import bruteForce, movingAverage

class TestUtil(unittest.TestCase):

    
    def testBruteForce(self):
        def isDifferent(l1, l2):
            for i in range(len(l1)):
                if l1[i] != l2[i]: return True
            return False
        elements = range(randint(1,7))
        count = randint(1,5)
        combList = [i for i in bruteForce(count, elements)]
        self.assertEqual(len(combList),len(elements)**count)
        for i in range(len(combList)):
            for j in range(len(combList)):
                if i != j:
                    self.assertTrue(isDifferent(combList[i], combList[j]))
    def testMovingAverage(self):
        #movingAverage([40, 30, 50, 46, 39, 44]) --> 40.0 42.0 45.0 43.0
        resList = [40.0, 40.0, 41.0, 41.8, 43, 44]
        mvAvg = movingAverage([40, 30, 50, 46, 39, 44])
        for res in resList:
            self.assertAlmostEqual(mvAvg.next(), res)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test']
    unittest.main()