import processes as proc
import translation
import molecules as mol


class Model:
    """
    Initializes the states and processes for the model and lets the processes update their corresponding states.
    """

    def __init__(self):
        self.states = {}
        self.processes = {}

        # initiate states
        self.ribosomes = {'Ribosomes': mol.Ribosome('Ribosomes', 'Ribosomes', 10)}
        self.mrnas = {'MRNA_{0}'.format(i): mol.MRNA(i, 'MRNA_{0}'.format(i), "UUUUUUUUUUAA") for i in range(50)}
        self.states.update(self.ribosomes)
        self.states.update(self.mrnas)

        # initiate processes
        trsl = translation.Translation(1, "Translation")
        trsl.set_states(self.mrnas.keys(), self.ribosomes.keys())
        self.processes = {"Translation": trsl}

    def step(self):
        """
        Do one update step for each process.

        """
        for p in self.processes:
            self.processes[p].update(self)

    def simulate(self, steps, log=True):
        """
        Simulate the model for some time.

        """
        for s in range(steps):
            self.step()
            if log:  # This could be an entry point for further logging
                # print count of each protein to the screen
                print('\r{}'.format([len(self.states[x]) for x in self.states.keys() if "Protein_" in x]), end='')


if __name__ == "__main__":
    c = Model()
    c.simulate(100, log=True)
