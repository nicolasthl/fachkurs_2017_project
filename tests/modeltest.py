import unittest
import modeldata
import molecules as mol
import model

class TestModel(unittest.TestCase):
    def setUp(self):
        self.m = model.Model()

    def test_timestep(self):
        '''
        calling step should increase the time
        @return:
        '''
        start = self.m.timestep
        self.m.step()
        self.assertGreater(self.m.timestep, start)

    def test_simulate_time(self):
        '''
        simulate should increase timestep by x
        @return:
        '''
        x = 100
        start = self.m.timestep
        self.m.simulate(x, log=False)
        self.assertEqual(start + x, self.m.timestep)


class TestData(unittest.TestCase):
    def test_data_mrna(self):
        '''
        calling get states with the mol.MRNA class should return ids with mRNA in them.
        @return:
        '''
        db = modeldata.ModelData()
        mrnas = db.get_states(mol.MRNA)
        for mrna in mrnas:
            self.assertRegex(mrna[1], "MRNA_\d+")


if __name__ == '__main__':
    unittest.main()