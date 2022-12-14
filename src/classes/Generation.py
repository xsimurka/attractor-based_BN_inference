from src.classes.BooleanNetwork import BooleanNetwork
from typing import Optional, List, Tuple
from scipy.special import softmax
import numpy as np
import pandas as pd
from random import choices
import src.classes.TargetBN as bn

State = Tuple[bool]


class Generation:
    """Class represents one generation of BNs of genetic algorithm

    Attributes:
        num_of_nets  number of networks in generation
        target_bn    object that holds the information about target BN
        networks     list of <num_of_nets> instances of <BooleanNetwork> that Generation consists of
        scores       fitness score for each network from <networks> in [0; 1]
                     describes how well the network fits the input data"""

    def __init__(self, num_of_nets: int, target_bn: bn.TargetBN,
                 nets: Optional[List[BooleanNetwork]] = None):
        self.num_of_nets = num_of_nets
        self.target_bn = target_bn
        if nets is None:
            self.networks = [BooleanNetwork(target_bn) for _ in range(num_of_nets)]
        else:
            self.networks = list(
                map(lambda n: n.deepcopy_(), nets))  # due to a possibility of picking the same 2 network instances
        self.scores = [0.0 for _ in range(num_of_nets)]

    @property
    def best(self):
        """Best score among all networks in generation"""

        return max(self.scores)

    def create_new_generation(self, num_of_mut_genes: int, num_of_mutations: int, elite_ratio: float) -> 'Generation':
        """Creates new generation of genetic algorithm by picking and mutating best BNs from the previous generation.

        :param num_of_mut_genes  number of genes to be mutated in each picked network
        :param num_of_mutations  number of mutations to be performed on each gene
        :param elite_ratio       real number from range [0;1] determining percentage of best fitting networks
                                 that are picked automatically to the next generation
        :return                  instance of new generation created"""

        scores = pd.Series(self.scores)
        num_of_elite = round(self.num_of_nets * elite_ratio)
        elite = scores.nlargest(num_of_elite).index.values.tolist()  # indices of best nets in generation
        weights = list(softmax(np.array(self.scores)))  # probability distribution of picking nets to new generation
        # by setting <k> argument as follows, the new generations will contain the same number of nets as previous one
        picked = choices(range(self.num_of_nets), weights=weights, k=self.num_of_nets - num_of_elite - 1)
        new_gen = Generation(self.num_of_nets, self.target_bn,
                             [self.networks[elite[0]]] + [self.networks[i] for i in elite] + [self.networks[i] for i in
                                                                                              picked])
        new_gen.mutate(num_of_mut_genes, num_of_mutations, 1)
        return new_gen

    def mutate(self, num_of_genes: int, num_of_mutations: int, begin: int) -> None:
        """Mutates given number of genes of each BN of given generation by given number of mutations

        :param num_of_genes      number of genes to be mutated
        :param num_of_mutations  number of mutations to be performed on each gene
        :param begin             index of first network that will be mutated"""

        for i in range(begin, len(self.networks), 1):  # skip the best network from the previous generation
            self.networks[i].mutate(num_of_genes, num_of_mutations)

    def compute_fitness(self):
        """Computes fitness of each BN of given generation comparing to target network described by the input data.
        Sets .scores attribute for each BN in generation"""

        for i in range(self.num_of_nets):
            total_score = self.networks[i].compute_fitness(-1)  # wild-type
            for j in sorted(self.target_bn.ko_sinks.keys()):  # knock-outs
                total_score += self.networks[i].compute_fitness(j, False)

            for j in sorted(self.target_bn.oe_sinks.keys()):  # over-expressions
                total_score += self.networks[i].compute_fitness(j, True)

            # final net's score is average score of all experiments
            self.scores[i] = total_score / (len(self.target_bn.ko_sinks) + len(self.target_bn.oe_sinks) + 1)
