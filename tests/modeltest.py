import unittest
import modeldata
import molecules as mol
import model

class TestModel(unittest.TestCase):
    def test_data_mrna(self):
        '''
        calling get states with the mol.MRNA class should return ids with mRNA in them.
        @return:
        '''
        db = modeldata.ModelData()
        mrnas = db.get_states(mol.MRNA)
        for mrna in mrnas:
            self.assertRegex(mrna[1], "MRNA_\d+")

    def test_timestep(self):
        '''
        calling step should increase the time
        @return:
        '''
        m = model.Model()
        start = m.timestep
        m.step()
        self.assertGreater(m.timestep, start)

if __name__ == '__main__':
    unittest.main()