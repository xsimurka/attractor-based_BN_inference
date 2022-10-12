import src.utils as utils
from src.detect import detect_steady_states, get_model_computational_structures, is_attractor_state
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
        target_sinks      steady-states observed from input data, which are targeted by the inference
        networks          list of <num_of_nets> instances of <BooleanNetwork> that Generation consists of
        scores            fitness score for each network from <networks> in [0; 1]
                          describes how well the network fits the input data"""

    def __init__(self, num_of_nets: int, target_bn_info: bn.BNInfo,
                 nets: Optional[List[BooleanNetwork]] = None):
        self.num_of_nets = num_of_nets
        self.target_bn_info = target_bn_info

        if nets is None:
            self.networks = [BooleanNetwork(target_bn_info) for _ in range(num_of_nets)]
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
        new_gen = Generation(self.num_of_nets, self.target_bn_info,
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

            for j in sorted(self.target_bn_info.ko_sinks.keys()):
                total_score += self.evaluate_fitness(i, j, False)

            for j in sorted(
                    self.target_bn_info.oe_sinks.keys()):  # iterates over perturbed gene indices of all experiments
                total_score += self.evaluate_fitness(i, j, True)

            # final net's score is average score of all experiments
            self.scores[i] = total_score / (len(self.target_bn_info.ko_sinks) + len(self.target_bn_info.oe_sinks) + 1)

    def evaluate_fitness(self, network_index: int, perturbed_gene_index: int,
                         perturbed_gene_state: Optional[bool] = None) -> float:
        """Evaluates fitness of given BN depending on its steady-state attractors comparing to steady-states from
        given attractor data.

        :param network_index          index of currently evaluated BN in actual generation
        :param perturbed_gene_index   index of perturbed gene of currently evaluated BN (-1 if wild-type)
        :param perturbed_gene_state   True if current experiment is over-expression of <perturbed_gene_index>
                                      False if current experiment is knock-out of <perturbed_gene_index>
                                      None if current experiment is wild-type
        :return                       real number from [0;1] that determines BN's fitness to input data"""

        isolated_variables = self.networks[network_index].get_isolated_variables(perturbed_gene_index)
        aeon_model_string = self.networks[network_index].to_aeon_string(perturbed_gene_index, isolated_variables,
                                                                        perturbed_gene_state)
        model, sag = get_model_computational_structures(aeon_model_string)
        target_sinks = self.get_target_sinks(perturbed_gene_index, perturbed_gene_state)

        if isolated_variables:  # reduce dimension of target sinks in order to match observed data dimension
            target_sinks = utils.reduce_attractors_dimension(target_sinks, isolated_variables)

        observed_sinks = detect_steady_states(sag)
        bpg = bg.BipartiteGraph(target_sinks, observed_sinks)
        cost, pairs = bpg.minimal_weighted_assignment()
        dimension = self.target_bn_info.num_of_vars - len(isolated_variables)  # number of variables in each state after reduction
        matching_variables = 0  # variable = one element of a state

        # if some states were matched, then total number of matching variables is equal to total number of variable
        # pairs minus cost (total number of variables that do not match), else it stays equal to 0 (initial value)
        if cost is not None:
            matching_variables += (min(len(target_sinks), len(observed_sinks)) * dimension) - cost

        if len(target_sinks) > len(observed_sinks):
            # some steady-states from data were not reached by actual model, which is problematic because such
            # model can not capture desired behaviour specified by the input data
            # total_variables is therefore, equal to number of variable pairs plus number of unpaired variables
            total_variables = len(target_sinks) * dimension
            unmatched_states = utils.get_unmatched_states(bpg, pairs, target_sinks, position=0)
            # check if missing steady-states are on some other type of attractor because there is possibility of
            # breaking such attractor to steady-state by tiny mutation of actual BN
            for state in unmatched_states:
                if is_attractor_state(model, sag, state):
                    matching_variables += dimension * 1 / 2  # penalty 1/2 for not being in single state attractor

        else:  # len(target_sinks) <= len(observed_sinks)
            # there is possibility that some steady-states were not caught while measuring which is more and more
            # probable by increasing number of genes where the number of steady-states grows exponentially
            # therefore, no penalty is given and total number of variables is only equal to number of variable pairs
            total_variables = len(target_sinks) * dimension

        return matching_variables / total_variables

    def get_target_sinks(self, perturbed_gene_index: int, perturbed_gene_state: Optional[bool] = None) -> List[State]:
        """Returns list of states of specific experiment from input data.

        :param perturbed_gene_index  index of perturbed gene
        :param perturbed_gene_state  True if over-expression is desired
                                     False if knock-out is desired
                                     None if wild type is desired (combined with -1 value of <perturbed_gene_index>
        :return                      list of states of specified experiment from the input data"""

        if perturbed_gene_index == -1:
            return self.target_bn_info.wt_sinks

        return (self.target_bn_info.oe_sinks if perturbed_gene_state else self.target_bn_info.ko_sinks)[
            perturbed_gene_index]
