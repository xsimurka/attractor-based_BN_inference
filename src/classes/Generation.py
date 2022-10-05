import src.utils as utils
from src.detect import detect_steady_states, get_symbolic_async_graph, is_attractor_state
from src.classes.BooleanNetwork import BooleanNetwork
from typing import Optional, List, Tuple
from scipy.special import softmax
import numpy as np
import pandas as pd
from random import choices
from copy import deepcopy
import src.classes.BNInfo as bn
import src.classes.BipartiteGraph as bg

State = Tuple[bool]


class Generation:
    """Class represents one generation of BNs of genetic algorithm

    Attributes:
        num_of_nets       number of networks in generation
        num_of_variables  number of gene variables in each network
        target_sinks      steady-states observed from input data, which are targeted by the inference
        networks          list of <num_of_nets> instances of <BooleanNetwork> that Generation consists of
        scores            fitness score for each network from <networks> in [0; 1]
                          describes how well the network fits the input data"""

    def __init__(self, num_of_nets: int, num_of_variables: int, target_sinks: bn.BNInfo,
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

            total_score += self.evaluate_fitness(i, -1)

            for j in sorted(self.target_sinks.ko_sinks.keys()):
                total_score += self.evaluate_fitness(i, j, False)

            for j in sorted(
                    self.target_sinks.oe_sinks.keys()):  # iterates over perturbed gene indices of all experiments
                total_score += self.evaluate_fitness(i, j, True)

            # final net's score is average score of all experiments
            self.scores[i] = total_score / (len(self.target_sinks.ko_sinks) + len(self.target_sinks.oe_sinks) + 1)

    def evaluate_fitness(self, network_index: int, perturbed_gene_index: int,
                         perturbed_gene_state: Optional[bool] = None) -> float:
        """Evaluates fitness of given BN depending on its steady-state attractors comparing to steady-states from
        given attractor data.

        :param network_index          biodivine_aeon.BooleanNetwork model of actual BN
        :param perturbed_gene_index   Symbolic Asynchronous Graph of actual model
        :param perturbed_gene_state   steady-states of particular experiment from data
        :return                       real number from [0;1] that determines BN's fitness
        TODO: with increasing number of nodes, number of steady-states grows exponentially and metrics becomes unusable
        """

        isolated_variables = self.networks[network_index].get_isolated_variables(perturbed_gene_index)
        # if WT then no gene is perturbed, if MT then looks at j-th variable of first (but basically any)
        # steady-state and derives type of perturbation (if 0 then it is KO, if 1 then it is OE)
        aeon_model_string = self.networks[network_index].to_aeon_string(perturbed_gene_index, isolated_variables,
                                                                        perturbed_gene_state)
        model, sag = get_symbolic_async_graph(aeon_model_string)
        target_sinks = self.get_target_sinks(perturbed_gene_index, perturbed_gene_state)
        if isolated_variables:
            target_sinks = utils.reduce_attractors_dimension(target_sinks, isolated_variables)
        observed_sinks = detect_steady_states(sag)
        bpg = bg.BipartiteGraph(target_sinks, observed_sinks)
        cost, pairs = bpg.minimal_weighted_assignment()
        dimension = len(target_sinks[0])
        # weight of one variable in one state is computed as:
        # number of overlapping state tuples + number of overhanging states, multiplied by dimension of reduced model
        # therefore, each assigned tuple of states "behaves as one" and has weight equal to their reduced dimension and
        # each overhanging state alone has weight equal to its reduced dimension, one variable thus, have weight equal to
        # 1 / total number of variables where matching tuple of variables act as one
        matching_variables = 0
        total_variables = (abs(len(target_sinks) - len(observed_sinks)) + min(len(target_sinks),
                                                                              len(observed_sinks))) * dimension

        # if some states were matched, then total number of matching variables is equal to total number of variables
        # minus cost (variables that do not match), else it stays equal to 0 (initial value)
        if cost is not None:
            matching_variables += ((min(len(target_sinks), len(observed_sinks)) * dimension) - cost)

        # try to observe how many sinks absent and on which side
        if len(target_sinks) > len(observed_sinks):
            # not ideal - some steady-states from data were not reached by given model,
            # check if missing steady-states are on some other type of attractor
            unmatched_states = utils.get_unmatched_states(bpg, pairs, 0, len(target_sinks))
            for state in unmatched_states:
                if is_attractor_state(model, sag, state):
                    matching_variables += dimension * 1 / 2  # penalty 1/3 for not being in single state

        elif len(target_sinks) < len(observed_sinks):
            # there is possibility that some steady-states were not caught while measuring, not a big problem if only few
            unmatched_states = utils.get_unmatched_states(bpg, pairs, 1, len(observed_sinks))
            matching_variables += len(unmatched_states) * dimension * 3 / 4  # penalty for not being in data

        # if target == observed then no correction is needed
        return matching_variables / total_variables

    def get_target_sinks(self, perturbed_gene_index: int, perturbed_gene_state: Optional[bool] = None) -> List[State]:
        """"""

        if perturbed_gene_index == -1:
            return self.target_sinks.wt_sinks

        return (self.target_sinks.oe_sinks if perturbed_gene_state else self.target_sinks.ko_sinks)[
            perturbed_gene_index]
