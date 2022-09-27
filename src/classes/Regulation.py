class Regulation:
    """Class represents structured regulation.
    Attributes:
        source  gene that regulates <target> gene
        target  gene that is regulated by <source> gene
        sign    type of regulation - positive (True) or negative (False)"""

    def __init__(self, source: int, target: int, sign: bool):
        self.source = source
        self.target = target
        self.sign = sign
