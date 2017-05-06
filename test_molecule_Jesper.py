import pytest
from molecules_Jesper import *


def test_valid_molecule_success():
    assert isinstance(Protein('ProteinA'), Protein)
    assert isinstance(Protein('ProteinB', 'ARND'), Protein)

    assert isinstance(Ribo('free'), Ribo)

    assert isinstance(MRNA('mRNA1'), MRNA)
    assert isinstance(MRNA('mRNA2', 'ACU'), MRNA)


def test_invalid_molecule_fails():
    with pytest.raises(ValueError):
        Protein('ProteinC', 'WX')

    with pytest.raises(ValueError):
        MRNA('mRNA3', 'TTT')


def test_equality():
    assert Protein('ProtA') == Protein('ProtA')
    assert Protein('ProtA') != Protein('ProtB')
    assert Protein('ProtA', 'ARND') != Protein('ProtA', 'DNRA')


def test_length():
    assert len(Protein('ProtA')) == 0
    assert len(Protein('ProtA', 'ARND')) == 4


def test_population_collection():
    ribo_coll = PopulationCollection(Ribo)

    with pytest.raises(ValueError):
        ribo_coll.add(Protein('Bla'))

    ribo_coll.add(Ribo('free'), 1000)
    ribo_coll.add(Ribo('bound'), 500)

    ribo_coll.pop('free', 100)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 500

    ribo_coll.pop('bound', 50)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 450


def test_particle_collection_basics():
    """Exactly same test as "test_population_collection", but now for
    a ParticleCollcetion. Asserts interface is identical."""
    ribo_coll = ParticleCollection(Ribo)

    with pytest.raises(ValueError):
        ribo_coll.add(Protein('Bla'))

    ribo_coll.add(Ribo('free'), 1000)
    ribo_coll.add(Ribo('bound'), 500)

    ribo_coll.pop('free', 100)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 500

    ribo_coll.pop('bound', 50)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 450


def test_particle_collection_pop():
    """The pop method for ParticleCollection should return the molecules."""
    prot_coll = ParticleCollection(Protein)
    prot_coll.add(Protein('ProtA'), 10)

    popped = prot_coll.pop('ProtA')
    assert len(popped) == 1
    assert prot_coll.count('ProtA') == 9


def test_particle_collection_count_match():
    """The count method for ParticleCollection allows matching on molecule
    properties such as polymer length."""
    prot_coll = ParticleCollection(Protein)
    prot_coll.add(Protein('ProtB', 'AGPRHS'), 15)

    popped = prot_coll.pop('ProtB', 5)
    assert len(popped) == 5
    assert prot_coll.count('ProtB') == 10

    for p in popped:
        p.add_monomer('S')
        prot_coll.add(p)

    assert prot_coll.count('ProtB') == 15
    assert prot_coll.count('ProtB', lambda p: len(p) == 6) == 10
    assert prot_coll.count('ProtB', lambda p: len(p) == 7) == 5


def test_particle_collection_pop_match():
    """The pop method for ParticleCollection allows matching on molecule
    properties such as polymer length."""
    prot_coll = ParticleCollection(Protein)
    prot_coll.add(Protein('ProtB', 'AGP'), 5)
    prot_coll.add(Protein('ProtB', 'AGPRHS'), 15)

    popped = prot_coll.pop('ProtB', number=0, matcher=lambda p: len(p) > 5)
    assert len(popped) == 15
    assert prot_coll.count('ProtB') == 5
