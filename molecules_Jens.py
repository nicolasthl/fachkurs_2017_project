class MoleculeContainer:
    def __init__(self, name, count):
        """
        
        Parameters
        ----------
        name : str
        count : int
        """
        self.name = name
        self.count = [self.name] * count

    def __repr__(self):
        return self.name


class Polymer:
    def __init__(self, name, id, sequence, valid_monomers):
        """
        Parameters
        ----------
        name : str
        id : int
        sequence : str
        valid_monomers : dict
        """
        self.id = id
        self.sequence = sequence
        self.name = name
        self.valid_monomers = valid_monomers#

    def __repr__(self):
        return self.name + '_' + str(self.id)

    def add_monomer(self, monomer):
        """
        Parameters
        ----------
        monomer : str
        """
        self.sequence += monomer

    def calc_mass(self):
        return sum(self.valid_monomers[monomer] for monomer in self.monomers)


class Ribo(MoleculeContainer):
    def __init__(self, name, count):
        MoleculeContainer.__init__(self, name, count)


class Protein(Polymer):
    amino_acids = {
        'A': 89.0929, 'R': 175.208, 'N': 132.118, 'D': 132.094, 'C': 121.158, 'Q': 146.144, 'E': 146.121,
        'G': 75.0664, 'H': 155.154, 'I': 131.172, 'L': 131.172, 'K': 147.195, 'M': 149.211, 'F': 165.189,
        'P': 115.13,  'S': 105.092, 'T': 119.119, 'W': 204.225, 'Y': 181.188, 'V': 117.146
    }  # comes from model data base

    def __init__(self, name, id, sequence=''):
        Polymer.__init__(self, name, id, sequence, self.amino_acids)


class mRNA(Polymer):
    nuc_acids = {'A': 1.0, 'U': 2.2, 'G': 2.1, 'C': 1.3}  # comes from model data base

    def __init__(self, name, id, sequence=''):
        Polymer.__init__(self, name, id, sequence, self.nuc_acids)
        self.sequence_triplet_binding = [0] * (len(self.sequence) // 3)


if __name__ == "__main__":
    ribo = Ribo('ribosome', 20)
    protein = Protein('Clb5', 1, "AYFG")
    mrna = mRNA('CLB5', 1, 'AGTTTCAAACCC')
    print(ribo)
    print(protein)
    print(mrna)