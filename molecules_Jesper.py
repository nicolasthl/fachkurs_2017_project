from collections import defaultdict
from copy import copy


class Molecule:
    def __init__(self, name):
        self.name = name


class Polymer:
    def __init__(self, valid_monomers, monomers):
        self.valid_monomers = valid_monomers
        self.monomers = ''
        for monomer in monomers:
            self.add_monomer(monomer)

    def __len__(self):
        return len(self.monomers)

    def add_monomer(self, monomer):
        if monomer not in self.valid_monomers:
            raise ValueError('Invalid monomer {}'.format(monomer))
        self.monomers += monomer

    def calc_mass(self):
        return sum(self.valid_monomers[monomer] for monomer in self.monomers)


class Ribo(Molecule):
    def __eq__(self, other):
        if not isinstance(other, Ribo):
            return NotImplemented
        else:
            return self.name == other.name


class Protein(Molecule, Polymer):
    amino_acids = {
        'A': 89.0929, 'R': 175.208, 'N': 132.118, 'D': 132.094, 'C': 121.158, 'Q': 146.144, 'E': 146.121,
        'G': 75.0664, 'H': 155.154, 'I': 131.172, 'L': 131.172, 'K': 147.195, 'M': 149.211, 'F': 165.189,
        'P': 115.13,  'S': 105.092, 'T': 119.119, 'W': 204.225, 'Y': 181.188, 'V': 117.146
    }

    def __init__(self, name, monomers=''):
        Molecule.__init__(self, name)
        Polymer.__init__(self, self.amino_acids, monomers)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            return NotImplemented
        else:
            return self.name == other.name and self.monomers == other.monomers


class MRNA(Molecule, Polymer):
    nuc_acids = {'A': 1.0, 'U': 2.2, 'G': 2.1, 'C': 1.3}

    def __init__(self, name, monomers=''):
        Molecule.__init__(self, name)
        Polymer.__init__(self, self.nuc_acids, monomers)

    def __eq__(self, other):
        if not isinstance(other, MRNA):
            return NotImplemented
        else:
            return self.name == other.name and self.monomers == other.monomers


class MoleculeCollection:
    def __init__(self, molecule_type):
        self.molecule_type = molecule_type
        self.molecules = None

    def add(self, molecule, number):
        """Adds number molecules to the collection."""
        if not isinstance(molecule, self.molecule_type):
            raise ValueError('Expected object of type {}, received of type {}'
                             .format(self.molecule_type, type(molecule)))

    def pop(self, name, number):
        """Removes number molecules named name from the collection. When number equals
        zero, all are removed."""
        pass

    def count(self, name):
        """Counts molecules named name."""
        pass


class PopulationCollection(MoleculeCollection):
    def __init__(self, molecule_type):
        super().__init__(molecule_type)
        self.molecules = defaultdict(int)

    def add(self, molecule, number=1):
        super().add(molecule, number)
        self.molecules[molecule.name] += number

    def pop(self, name, number=1):
        if number == 0:
            self.molecules[name] = 0
        else:
            self.molecules[name] -= number

    def count(self, name):
        return self.molecules[name]


class ParticleCollection(MoleculeCollection):
    def __init__(self, molecule_type):
        super().__init__(molecule_type)
        self.molecules = defaultdict(list)

    def add(self, molecule, number=1):
        super().add(molecule, number)
        for _ in range(number):
            self.molecules[molecule.name].append(copy(molecule))

    def pop(self, name, number=1, matcher=lambda x: True):
        """To select for molecules with certain properties, such as polymer length,
        provide a matcher function that evaluates to True for matching molecules."""
        result = []
        non_matching = []
        finished = False

        while self.molecules[name] and not finished:
            molecule = self.molecules[name].pop()
            if matcher(molecule):
                result.append(molecule)
            else:
                non_matching.append(molecule)

            finished = False if number == 0 else len(result) == number

        assert number == 0 or len(result) == number
        self.molecules[name] += non_matching

        return result

    def count(self, name, matcher=lambda x: True):
        return len([x for x in self.molecules[name] if matcher(x)])
