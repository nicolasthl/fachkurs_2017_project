import processes
from molecules import Ribo, Protein, MRNA, PopulationCollection, ParticleCollection


class Translation(processes.Process):
    """
    Translation is instantiated in the Cell to produce proteins.

    Defines Translation process. It iterates over all ribosomes and decides what
    they should do. They either bind to mRNA or elongate/terminate a protein if
    they are already bound.

    """

    def __init__(self, name, model):
        # call the constructor of the base class (processes.Process in this case)
        super().__init__(name, model)

    def __str__(self):
        # return string output of translation process 
        # todo: each process class should define this
        return "Translation process for mRNAs: {}".format(list(self.model.states['mRNA']))

    def update(self):
        """
        Update all mrnas and translate proteins.
        """
        for mrna_id in self.model.states[MRNA].molecules:
            prot = None
            for mrna in self.model.states[MRNA].molecules[mrna_id]:
                self.initiate(mrna)
                self.elongate(mrna)

    def initiate(self, mrna):
        """
        Try to bind to a given MRNA. Binding probability corresponds to the ribosome count.

        @type mrna: MRNA
        """
        pass

    def elongate(self, mrna):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.
        """
        pass

    def terminate(self, mrna):
        """
        Splits the ribosome/MRNA complex and returns a protein.
        """
        pass
