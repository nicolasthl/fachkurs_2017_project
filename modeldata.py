import pandas
import sqlite3
import os.path
import random as rnd
import molecules as mol

class ModelData:
    """
    class to process data sources for usage in the model
    """
    def __init__(self):
        pass

    def get_states(self, molecule_class):
        """
        retrieves the information required to construct the different model molecules
        @param molecule_class: BioMolecule class
        @return: list
        """
        if molecule_class == mol.MRNA:
            alphabet = "AUGC"
            sequence = ''.join([rnd.choice(alphabet) for i in range(rnd.randint(50,500))])
            mrnas = []
            for i in range(50):
                 mrnas.append(('MRNA_{0}'.format(i), sequence))
            return mrnas
