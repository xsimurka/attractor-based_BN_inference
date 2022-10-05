from random import random


class Regulation:
    """Class represents structured regulation.
    Attributes:
        source  gene that regulates <target> gene
        target  gene that is regulated by <source> gene
        sign    type of regulation - positive (True) or negative (False)"""

    def __init__(self, gene1: int, gene2: int, sign: bool, directed: bool):
        if directed or random() < 0.5:
            self.source = gene1
            self.target = gene2
        else:
            self.source = gene2
            self.target = gene1

        self.sign = sign
