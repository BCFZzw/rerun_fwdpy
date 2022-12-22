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
        """
        genotype = np.array([[i, i] for i in range(0, 51)])
        pos = np.array([i for i in range(0, 51)])

        maxLen = 50
        boundaries = 25

        segment = ws.segmentation(genotype, pos, maxLen, boundaries)
        print(segment)


                  



if __name__ == '__main__':
    unittest.main()