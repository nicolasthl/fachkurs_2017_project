import modeldata
import molecules
import translation
import processes


class Output:
    """
    class for handling the simulation results of the different species types
    """

    def __init__(self, model):
        self.model = model
        self.timecourses = {state: SimulationResult(model.states[state]) for state in model.states}

    def add_timepoint(self, species):
        """
        add a simulation time point for one species
        @param species: mol.BioMolecule
        @return: None
        """
        if isinstance(self.model.states[species], molecules.BioMoleculeSet):
            pass  # TODO: implement a useful method BioMoleculeSet
        elif isinstance(self.model.states[species], molecules.BioMoleculeCount):
            self.timecourses[species].add_timepoint(self.model.states[species].count, self.model.timestep)


class SimulationResult:
    """
    handles and stores a simulation result for one species
    """

    def __init__(self, species):
        self.name = species.name
        self.value = []
        self.time = []

    def add_timepoint(self, time, value):
        self.value.append(value)
        self.time.append(time)


class Model:
    """
    Initializes the states and processes for the model and lets the processes update their corresponding states.
    """

    def __init__(self):
        # all selfs should be initialized in the constructor
        self.states = {}
        self.processes = {}
        self.timestep = 0
        self.db = modeldata.ModelData()

        self.__initialize_ribosomes()
        self.__initialize_mRNA()
        self.__initialize_processes()

        self.results = Output(self)

    def add_new_state(self, molecule_container):
        assert isinstance(molecule_container, molecules.BioMoleculeContainer)
        self.states[molecule_container.name] = molecule_container

    def add_new_process(self, process):
        assert isinstance(process, processes.Process)
        self.processes[process.name] = process

    def __initialize_ribosomes(self):
        self.add_new_state(molecules.Ribosome('Ribosomes', 10))

    def __initialize_mRNA(self):
        self.add_new_state(molecules.BioMoleculeSet('mRNA'))
        for mid, sequence in self.db.get_states(molecules.MRNA):
            self.states['mRNA'][mid] = molecules.MRNA(mid, sequence)

    def __initialize_processes(self):
        trsl = translation.Translation("Translation", self)
        self.add_new_process(trsl)

    def step(self):
        """
        Do one update step for each process.

        """
        for p in self.processes:
            self.processes[p].update()

        for state in self.states:
            self.results.add_timepoint(state)

        self.timestep += 1

    def simulate(self, steps, log=True):
        """
        Simulate the model for some time.

        """
        for s in range(steps):
            self.step()
            if log:  # This could be an entry point for further logging
                # print count of each protein to the screen
                print('{}'.format([(len(self.states[x]), x) for x in self.states.keys() if "Protein" in x]))


if __name__ == "__main__":
    c = Model()
    c.simulate(100, log=True)
