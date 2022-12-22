import unittest
import windowStats as ws
import pandas as pd
import allel
import numpy as np
import os



class TestParsingMethods(unittest.TestCase):
    
    def test_segmentation(self):
        """
        The previous version will miss the starting region in the middle 
        The beginning and the end is good
        but it shouldn't double the score
        fix it and run again
        """
        genotype = np.array([[i, i] for i in range(0, 51)])
        pos = np.array([i for i in range(0, 51)])

        maxLen = 50
        boundaries = 25

        segment = ws.segmentation(genotype, pos, maxLen, boundaries)
        print(segment)

    def test_shapeTransform(self):
        """
        The previous version uses np.sum, which works for simulated data 
        since there are no missing data, change to allel function instead
        """
        genotype = np.random.randint(2, size = [3, 3, 2])
        genotype = allel.GenotypeArray(genotype, dtype='i1')

        print(genotype)
        genotype012 = ws.shapeTransform(genotype)
        print(genotype012)


                  



if __name__ == '__main__':
    unittest.main()