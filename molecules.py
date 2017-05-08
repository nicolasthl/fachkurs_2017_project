from copy import copy


class Molecule:
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.name == other.name


class Polymer(Molecule):
    def __init__(self, name, sequence, valid_monomers):
        super().__init__(name)
        self.valid_monomers = valid_monomers
        self.sequence = ''
        self.bindings = []
        for monomer in sequence:
            self.add_monomer(monomer)

    def __eq__(self, other):
        return super().__eq__(other) and self.sequence == other.sequence

    def __len__(self):
        return len(self.sequence)

    def add_monomer(self, monomer):
        if monomer not in self.valid_monomers:
            raise ValueError('Invalid monomer {}'.format(monomer))
        self.sequence += monomer
        self.bindings.append(None)

    def calc_mass(self):
        return sum(self.valid_monomers[monomer] for monomer in self.sequence)


class Ribo(Molecule):
    pass


class Protein(Polymer):
    amino_acids = {
        'A': 89.0929, 'R': 175.208, 'N': 132.118, 'D': 132.094, 'C': 121.158, 'Q': 146.144, 'E': 146.121,
        'G': 75.0664, 'H': 155.154, 'I': 131.172, 'L': 131.172, 'K': 147.195, 'M': 149.211, 'F': 165.189,
        'P': 115.13,  'S': 105.092, 'T': 119.119, 'W': 204.225, 'Y': 181.188, 'V': 117.146
    }

    def __init__(self, name, sequence=''):
        super().__init__(name, sequence, self.amino_acids)


class MRNA(Polymer):
    nuc_acids = {'A': 1.0, 'U': 2.2, 'G': 2.1, 'C': 1.3}

    def __init__(self, name, sequence=''):
        super().__init__(name, sequence, self.nuc_acids)


class MoleculeCollection:
    def __init__(self, molecule_type):
        self.molecule_type = molecule_type
        self.molecules = None

    def add(self, molecule):
        if not isinstance(molecule, self.molecule_type):
            raise ValueError('Expected object of type {}, received of type {}'
                             .format(self.molecule_type, type(molecule)))

    def take(self, name, number):
        pass

    def count(self, name):
        pass

    def populate(self, name, number):
        for _ in range(number):
            self.add(self.molecule_type(name))


class PopulationCollection(MoleculeCollection):
    def __init__(self, molecule_type):
        super().__init__(molecule_type)
        self.molecules = dict()

    def add(self, molecule):
        super().add(molecule)
        if molecule.name in self.molecules.keys():
            self.molecules[molecule.name] += 1
        else:
            self.molecules[molecule.name] = 1

    def take(self, name, number=1):
        assert self.molecules[name] >= number
        self.molecules[name] -= number

    def count(self, name=None):
        if name is None:
            return sum([self.molecules[name] for name in self.molecules])
        return self.molecules[name]


class ParticleCollection(MoleculeCollection):
    def __init__(self, molecule_type):
        super().__init__(molecule_type)
        self.molecules = dict()

    def add(self, molecule):
        super().add(molecule)
        if molecule.name in self.molecules.keys():
            self.molecules[molecule.name].append(copy(molecule))
        else:
            self.molecules[molecule.name] = [copy(molecule)]

    def take(self, name, number=1):
        result = []
        while self.molecules[name]:
            molecule = self.molecules[name].pop()
            result.append(molecule)
            if number == len(result):
                break

        assert number == 0 or number == len(result)

        return result

    def count(self, name=None):
        if name is None:
            return sum([len([x for x in self.molecules[molname]]) for molname in self.molecules])
        return len([x for x in self.molecules[name]])

    def get_molecules(self, name=None):
        if not name:
            return [molecule for molecules in self.molecules.values() for molecule in molecules]
        else:
            return self.molecules[name]
