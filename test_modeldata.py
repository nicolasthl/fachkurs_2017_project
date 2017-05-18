import pytest
from database import *
from molecules import MRNA


def test_mrna_sequences():
    md = ModelData()
    mrnas = md.get_states(MRNA)
    for mrna_id, sequence in mrnas:
        for char in sequence:
            assert char in 'GUAC'

if __name__ == '__main__':
    pytest.main(['test_modeldata.py'])