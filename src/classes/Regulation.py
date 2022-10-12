class Regulation:
    """Class represents structured regulation.
    Attributes:
        source    gene that regulates <target> gene
        target    gene that is regulated by <source> gene
        sign      type of regulation - positive (True) or negative (False)"""

    def __init__(self, source: int, target: int, sign: bool):
        self.source = source
        self.target = target
        self.sign = sign

    def __repr__(self):
        return 'v_{} -{} v_{}'.format(self.source, ">" if self.sign else "|", self.target)

    def __eq__(self, other: 'Regulation'):
        return self.source == other.source and self.target == other.target and self.sign == other.sign

    def __hash__(self):
        return hash(repr(self))
