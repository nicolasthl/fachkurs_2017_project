class Process:
    """
    Parent for all cellular processes.
    """

    def __init__(self, name, model):
        self.name = name
        self.model = model

    def update(self):
        """
        Has to be implemented by child class.
        """
        pass
