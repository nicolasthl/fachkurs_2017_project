import pytest
from molecules import *


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

    ribo_coll.populate('free', 1000)
    ribo_coll.populate('bound', 500)

    ribo_coll.take('free', 100)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 500

    ribo_coll.take('bound', 50)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 450


def test_particle_collection_basics():
    """Exactly same test as "test_population_collection", but now for
    a ParticleCollcetion. Asserts interface is identical."""
    ribo_coll = ParticleCollection(Ribo)

    with pytest.raises(ValueError):
        ribo_coll.add(Protein('Bla'))

    ribo_coll.populate('free', 1000)
    ribo_coll.populate('bound', 500)

    ribo_coll.take('free', 100)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 500

    ribo_coll.take('bound', 50)
    assert ribo_coll.count('free') == 900
    assert ribo_coll.count('bound') == 450


def test_particle_collection_take():
    """The pop method for ParticleCollection should return the molecules."""
    prot_coll = ParticleCollection(Protein)
    prot_coll.populate('ProtA', 10)

    popped = prot_coll.take('ProtA')
    assert popped == [Protein('ProtA')]
    assert prot_coll.count('ProtA') == 9


def test_protein_iterate():
    p = Protein('Bla', 'ADRA')
    assert len([monomer for monomer in p.sequence]) == 4


def test_particle_iterate():
    prot_coll = ParticleCollection(Protein)
    prot_coll.populate('ProtA', 10)
    prot_coll.populate('ProtB', 20)

    assert len([p for p in prot_coll.get_molecules('ProtA')]) == 10
    assert len([p for p in prot_coll.get_molecules('ProtB')]) == 20
    assert len([p for p in prot_coll.get_molecules()]) == 30

if __name__ == '__main__':
    pytest.main(['test_molecules.py'])