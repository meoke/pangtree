class NoConsensusFound(Exception):
    pass


class NoTresholdFound(Exception):
    pass


class StopExceeded(Exception):
    pass


class CloseProgram(Exception):
    def __init__(self, message):
        self.message = message