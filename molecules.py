class BioMolecule:
    def __init__(self, mid, mass):
        self.mid = mid
        self.mass = mass


class BioMoleculeContainer:
    def __init__(self, name):
        self.name = name


class BioMoleculeCount(BioMoleculeContainer):
    def __init__(self, name, count=0):
        super().__init__(name) # call __init__ of BioMoleculeContainer
        self.count = count


class Ribosome(BioMoleculeCount):
    """
    A ribosome can bind MRNA and translate it. After translation is
    finished it produces a protein.

    During initiation the ribosome checks if a given MRNA is bound
    by another ribosome and binds only if position 0 is empty.

    Elongation checks if the next codon is unbound and only elongates
    if the ribosome can move on. If the ribosome encounters a stop
    codon ("*") translation terminates. The MRNA is detached from the
    ribosome and the finished protein is returned.
    """

    def __init__(self, name, count=0):
        super().__init__(name, count)


class Polymerase(BioMoleculeCount):
    """
    A polymerase that can elongate nucleotide molecules. This could be used to derive special
    RNA and DNA polymerases.
    """
    pass


class RNAPolymeraseII(Polymerase):
    """
    A polymerase that generates mRNAs from DNA sequences.
    """
    pass


class BioMoleculeSet(BioMoleculeContainer):
    def __init__(self, name, biomlist=[]):
        super().__init__(name)
        self.biomolecule_dict = {}
        for biom in biomlist:
            self.biomolecule_dict[biom.mid] = biom

    def __getitem__(self, key):
        return self.biomolecule_dict[key]

    def __setitem__(self, key, biom):
        if not isinstance(biom, BioMolecule):
            raise Exception('BioMoleculeList can only contain Biomolecules')
        assert biom.mid == key
        self.biomolecule_dict[key] = biom

    def __iter__(self):
        for key in self.biomolecule_dict.keys():
            yield key

    @property
    def count(self):
        return len(self.biomolecule_dict)


class Polymer(BioMolecule):
    """
    A polymer molecule that has a sequence attribute which is
    accessible via indexing the object.

    @type mid: str
    @type sequence: list or str
    @type mass: float
    """

    def __init__(self, mid, sequence, mass=0):
        super().__init__(mid, mass)
        self.sequence = sequence

    def __getitem__(self, key):
        """
        implement the access operator []
        @param key: int or slice
        @return: sequence element
        """
        return self.sequence[key]


class MRNA(Polymer):
    def __init__(self, mid, sequence, mass=0):
        super().__init__(mid, sequence, mass)
        self.sequence_triplet_binding = [0] * (len(sequence) // 3)
        self.calculate_mass()

    def calculate_mass(self):
        self.mass = 0
        NA_mass = {'A': 1.0, 'U': 2.2, 'G': 2.1, 'C': 1.3}
        for na in self.sequence:
            self.mass += NA_mass[na]


class Protein(Polymer):
    """
    Protein with Polymer features and mass calculation. A global class
    attribute counts the number of proteins that have been instantiated.

    A protein can be elongated using the "+" operator:

    >> protein = Protein(1, "prot", "MVFT")
    >> protein + "A"
    >> protein.sequence
    MVFTA
    """
    number_of_proteins = 0

    def __init__(self, mid, sequence, mass=0):
        super().__init__(mid, sequence, mass)

    def __iadd__(self, AS):
        self.sequence = self.sequence + AS
        return self

    def calculate_mass(self):
        AA_mass = dict(A=89.0929, R=175.208, N=132.118, D=132.094, C=121.158, Q=146.144, E=146.121, G=75.0664,
                       H=155.154, I=131.172, L=131.172, K=147.195, M=149.211, F=165.189, P=115.13, S=105.092, T=119.119,
                       W=204.225, Y=181.188, V=117.146)
        for aa in self.sequence:
            self.mass += AA_mass[aa]


