import pytest
from model import *


def test_timestep():
    """
    calling step should increase the time
    """
    m = Model()
    start = m.timestep
    m.step()
    assert m.timestep > start


def test_simulate_time():
    m = Model()
    x = 100
    start = m.timestep
    m.simulate(x, log=False)
    assert (start + x) == m.timestep

if __name__ == '__main__':
    pytest.main(['test_model.py'])