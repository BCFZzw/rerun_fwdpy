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

    def test_pandasParsing(self):
        """
        still using pandas as the new pipeline, if snp occurs at the last position, it will ignore it
        """
        
        pos = np.array([i for i in range(0, 10)])
        maxLen = 100
        binSize = 25

        df, indiceList = ws.pandasParsing(pos, maxLen, binSize)
        print(len(indiceList))


    def test_correctSegment(self):
        pos = np.array([i for i in range(0, 10)])
        genotype = np.random.randint(2, size = [10, 3, 2])
        genotype = allel.GenotypeArray(genotype, dtype='i1')
        genotype012 = ws.shapeTransform(genotype)
        maxLen = 100
        binSize = 25

        df, indiceList = ws.pandasParsing(pos, maxLen, binSize)
        
        index = 0
        
        segment = ws.correctSegment(genotype012, indiceList, index)
        print(segment)

        index = 1
        segment = ws.correctSegment(genotype012, indiceList, index)
        print(segment)

        index = 2
        segment = ws.correctSegment(genotype012, indiceList, index)
        print(segment)
        print(np.shape(segment))




                  



if __name__ == '__main__':
    unittest.main()