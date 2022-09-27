import src.utils as utils
from src.detect import detect_steady_states, get_symbolic_async_graph
from src.classes.BooleanNetwork import BooleanNetwork
from typing import Optional, List
from scipy.special import softmax
import numpy as np
import pandas as pd
from random import choices
from copy import deepcopy
from src.parse_input import BNInfo


class Generation:
    """Class represents one generation of BNs of genetic algorithm

    Attributes:
        num_of_nets       number of networks in generation
        num_of_variables  number of gene variables in each network
        target_sinks      steady-states observed from input data, which are targeted by the inference
        networks          list of <num_of_nets> instances of <BooleanNetwork> that Generation consists of
        scores            fitness score for each network from <networks> in [0; 1]
                          describes how well the network fits the input data"""

    def __init__(self, num_of_nets: int, num_of_variables: int, target_sinks: BNInfo,
                 nets: Optional[List[BooleanNetwork]] = None):
        self.num_of_nets = num_of_nets
        self.num_of_variables = num_of_variables
        self.target_sinks = target_sinks
        if nets is None:
            self.networks = [BooleanNetwork(num_of_variables) for _ in range(num_of_nets)]
        else:
            self.networks = list(map(deepcopy, nets))  # due to a possibility of picking the same 2 network instances
        self.scores = [0.0 for _ in range(num_of_nets)]

    @property
    def best(self):
        """Best score among all networks in generation"""

        return max(self.scores)

    def create_new_generation(self, num_of_mut_genes: int, num_of_mutations: int, best_ratio: float) -> 'Generation':
        """Creates new generation of genetic algorithm by picking and mutating best BNs from the previous generation.

        :param num_of_mut_genes  number of genes to be mutated in each picked network
        :param num_of_mutations  number of mutations to be performed on each gene
        :param best_ratio        real number from range [0;1] determining percentage of best fitting networks
                                 that are picked automatically to the next generation
        :return                  instance of new generation created"""

        num_of_best = round(self.num_of_nets * best_ratio)
        best = pd.Series(self.scores).nlargest(num_of_best).index.values.tolist()  # indices of best nets in generation
        weights = list(softmax(np.array(self.scores)))  # probability distribution of picking nets to new generation
        # by setting <k> argument as follows, the new generations will contain the same number of nets as previous one
        picked = choices(range(self.num_of_nets), weights=weights, k=self.num_of_nets - num_of_best)
        new_gen = Generation(self.num_of_nets, self.num_of_variables, self.target_sinks,
                             [self.networks[i] for i in best] + [self.networks[i] for i in picked])
        new_gen.mutate(num_of_mut_genes, num_of_mutations)
        return new_gen

    def mutate(self, num_of_genes: int, num_of_mutations: int) -> None:
        """Mutates given number of genes of each BN of given generation by given number of mutations

        :param num_of_genes      number of genes to be mutated
        :param num_of_mutations  number of mutations to be performed on each gene"""

        for net in self.networks:
            net.mutate(num_of_genes, num_of_mutations)

    def compute_fitness(self):
        """Computes fitness of each BN of given generation comparing to target network described by the input data.
        Sets .scores attribute for each BN in generation"""

        for i in range(self.num_of_nets):
            total_score = 0

            isolated_variables = self.networks[i].get_isolated_variables(-1)
            aeon_model_string = self.networks[i].to_aeon_string(-1, isolated_variables)
            model, sag = get_symbolic_async_graph(aeon_model_string)
            observed_sinks = detect_steady_states(sag)
            total_score += utils.evaluate_fitness(model, sag,
                                                  utils.reduce_attractors_dimension(self.target_sinks.wt_sinks,
                                                                                    isolated_variables),
                                                  observed_sinks)

            for j in sorted(self.target_sinks.ko_sinks.keys()):
                isolated_variables = self.networks[i].get_isolated_variables(j)
                aeon_model_string = self.networks[i].to_aeon_string(j, isolated_variables, False)
                model, sag = get_symbolic_async_graph(aeon_model_string)
                observed_sinks = detect_steady_states(sag)
                total_score += utils.evaluate_fitness(model, sag,
                                                      utils.reduce_attractors_dimension(self.target_sinks.ko_sinks[j],
                                                                                        isolated_variables),
                                                      observed_sinks)

            for j in sorted(
                    self.target_sinks.oe_sinks.keys()):  # iterates over perturbed gene indices of all experiments
                isolated_variables = self.networks[i].get_isolated_variables(j)
                # if WT then no gene is perturbed, if MT then looks at j-th variable of first (but basically any)
                # steady-state and derives type of perturbation (if 0 then it is KO, if 1 then it is OE)
                aeon_model_string = self.networks[i].to_aeon_string(j, isolated_variables, True)
                model, sag = get_symbolic_async_graph(aeon_model_string)
                observed_sinks = detect_steady_states(sag)
                total_score += utils.evaluate_fitness(model, sag,
                                                      utils.reduce_attractors_dimension(self.target_sinks.oe_sinks[j],
                                                                                        isolated_variables),
                                                      observed_sinks)
            # final net's score is average score of all experiments
            self.scores[i] = total_score / (len(self.target_sinks.ko_sinks) + len(self.target_sinks.oe_sinks) + 1)
