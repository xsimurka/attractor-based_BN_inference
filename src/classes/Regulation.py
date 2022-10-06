from random import random


class Regulation:
    """Class represents structured regulation.
    Attributes:
        source    gene that regulates <target> gene
        target    gene that is regulated by <source> gene
        sign      type of regulation - positive (True) or negative (False)
        directed  True if direction of regulation is fixed, otherwise False"""

    def __init__(self, gene1: int, gene2: int, sign: bool, directed: bool):
        if directed or random() < 0.5:
            self.source = gene1
            self.target = gene2
        else:  # if not directed and random value >= 0.5, then swap source and target gene
            self.source = gene2
            self.target = gene1

        self.sign = sign
