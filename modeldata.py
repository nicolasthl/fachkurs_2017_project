"""

"""
import pandas
import sqlite3
import os.path
import random as rnd
import molecules as mol

class ModelData:
    def __init__(self):
        pass

    def get_states(self, molecule_class):
        alphabet = "AUGC"
        sequence = ''.join([rnd.choice(alphabet) for i in range(50)])
        if molecule_class == mol.MRNA:
            mrnas = []
            for i in range(50):
                 mrnas.append((i, 'MRNA_{0}'.format(i), sequence))
            return mrnas
