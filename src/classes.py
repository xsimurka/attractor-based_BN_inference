from detect import detect_steady_states, get_symbolic_async_graph
from random import randint, choices, choice
from biodivine_aeon import SymbolicAsyncGraph
from typing import Set, Tuple, Optional, List, Dict, Any
from scipy.special import softmax
import numpy as np
import pandas as pd
from copy import deepcopy
from operator import itemgetter
from attractor_for_state import is_attractor_state
from scipy.optimize import linear_sum_assignment


State = Tuple[bool]


class Regulation:
    def __init__(self, source: int, target: int, sign: bool):
        self.source = source
        self.target = target
        self.sign = sign


class UpdateFunction:
    """Class represents boolean update rule in the form of Nested Canalyzing Function (NCF)

    Attributes:
        gene        id of regulated gene
        max_arity   total number of gene variables in BN i.e., maximum arity of the function
        canalyzing  canalyzing values
        canalyzed   canalyzed values
        indices     ids of genes that regulates <gene>
        fixed       set of fixed regulators of actual gene"""

    def __init__(self, gene: int, max_arity):
        self.gene = gene
        self.max_arity = max_arity
        self.canalyzing = []
        self.canalyzed = []
        self.indices = []
        self.fixed = set()

    @property
    def non_fixed(self) -> int:
        """Number of non fixed regulators"""
        return len(self.indices) - len(self.fixed)

    @property
    def arity(self) -> int:
        """Actual number of regulators"""
        return len(self.indices)


    def to_aeon_string(self) -> str:
        """Returns string representation of given NCF in .aeon format

        :return aeon string representation"""

        if self.arity == 0:  # an empty NCF formula can not occur in AEON string (results in an error)
            return str()

        # the most nested regulator will not have brackets around
        function_string = "{}v_{}".format("" if self.canalyzing[-1] == self.canalyzed[-1] else "!", self.indices[-1])

        for i in reversed(range(self.arity - 1)):  # without the last regulator already
            logic_operator = "&" if self.canalyzed[i] == 0 else "|"
            parity = "" if self.canalyzing[i] == self.canalyzed[i] else "!"
            function_string = "({}v_{} {} ".format(parity, self.indices[i], logic_operator) + function_string + ")"

        return "$v_{}:{}".format(self.gene, function_string)


    def add_gene(self, gene: int, fixed: bool, position: Optional[int] = None) -> None:
        """Adds new gene index to list of regulators of given NCF on given position. For time complexity reason it is
        differentiated between appending on the last position and inserting elsewhere.

        :param gene      id of added gene
        :param fixed     True if regulator can not be removed later during the inference, otherwise False
        :param position  specifies the position in NCF rule where to be inserted
                         if None, it is appended to the last position in NCF rule"""

        if position is None:
            self.indices.append(gene)
        else:
            self.indices.insert(position, gene)
        if fixed:
            self.fixed.add(gene)


    def add_rule(self, canalyzing: int, canalyzed: int, position: Optional[int] = None) -> None:
        """Adds given canalyzing and canalyzed value to given NCF rule on given position. For time complexity reason
        it is differentiated between appending on the last position and inserting elsewhere.

        :param canalyzing  canalyzing value
        :param canalyzed   canalyzed value
        :param position    specifies the position in NCF rule where to be inserted
                           if None, it is appended to the last position in NCF rule"""

        if position is None:
            self.canalyzing.append(canalyzing)
            self.canalyzed.append(canalyzed)
        else:
            self.canalyzing.insert(position, canalyzing)
            self.canalyzed.insert(position, canalyzed)


    def evaluate(self, interpretation: Tuple[int]) -> bool:
        """Returns the valuation of given NCF rule of given interpretation
        If variable from interpretation matches its corresponding canalyzing value, function returns corresponding
        canalyzed value. If no variable matches canalyzing value, functions evaluates itself to negation of the last
        canalyzed value.

        :param interpretation  Boolean vector specifying interpretation of atomic variables
        :return  Boolean valuation of given NCF using given interpretation"""

        assert (len(interpretation) == len(self.canalyzing) == len(self.canalyzed))
        for i in range(len(interpretation)):
            if interpretation[i] == self.canalyzing[i]:
                return bool(self.canalyzed[i])

        return not bool(self.canalyzed[-1])


    def select_possible_mutations(self, gene: int, total_num_of_regulations: int) -> List[str]:
        """Selects all reasonable mutations for particular selected gene according to its occurrence in NCF rule.
        Leaves possibility that no mutation is available - if all regulations of picked gene are fixed and model is
        already too dense.

        :param gene                      id of regulator to be checked
        :param total_num_of_regulations  the number of all current regulations in the network
        :return  list of possible mutations that are reasonable to be done on given gene"""

        mutations = []
        if gene in self.indices:
            if self.arity >= 2:
                mutations.append("c_and_c_values_swapping")

            if gene not in self.fixed:
                mutations.append("canalyzing_value_reversion")
                mutations.append("canalyzed_value_reversion")
                mutations.append("c_and_c_values_removal")

            if self.get_c_and_c_values(gene) in {(0, 1), (1, 0)}:
                mutations.append("c_and_c_values_reversion")

        elif total_num_of_regulations < (self.max_arity ** 2) / 2:  # restricted number of regulations in network
            mutations.append("c_and_c_values_insertion")

        return mutations


    def mutate(self, num_of_genes: int, num_of_mutations: int, total_num_of_regulations: int) -> None:
        """Selects randomly given amount of regulators to be mutated, select possible mutations of each
        regulators, performs chosen mutations

        :param num_of_genes              total number of genes in BN
        :param total_num_of_regulations  total number of regulations in the whole network
        :param num_of_mutations          number of mutations to be done"""

        regulators_to_mutate = choices(range(num_of_genes), k=num_of_mutations)

        for regulator in regulators_to_mutate:
            mutations = self.select_possible_mutations(regulator, total_num_of_regulations)
            if not mutations:
                continue
            getattr(self, choice(mutations))(regulator)


    def get_c_and_c_values(self, gene: int) -> Optional[Tuple[int, int]]:
        """Returns canalyzing and canalyzed value of particular gene id if occurred in given NCF rule

        :param gene  id of gene
        :return tuple with canalyzing and canalyzed value if gene occurs in given NCF rule, otherwise None"""

        try:
            i = self.indices.index(gene)
        except ValueError:
            return None
        return self.canalyzing[i], self.canalyzed[i]


    def get_nth_non_fixed_gene_index(self, n: int) -> int:
        """Returns id of nth non fixed regulator of given NCF rule

        :param n
        :return index of nth non fixed regulator, -1 if less than n regulators are non fixed"""

        result = -1
        if n > self.non_fixed:
            return result

        while n > 0:
            result += 1
            if self.indices[result] not in self.fixed:
                n -= 1

        return result


    def canalyzing_value_reversion(self, gene: int) -> None:
        """Reverts canalyzing value of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.indices.index(gene)
        self.canalyzing[i] = 1 - self.canalyzing[i]


    def canalyzed_value_reversion(self, gene: int) -> None:
        """Reverts canalyzed value of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.indices.index(gene)
        self.canalyzed[i] = 1 - self.canalyzed[i]


    def c_and_c_values_reversion(self, gene: int) -> None:
        """Reverts canalyzing and canalyzed values of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.indices.index(gene)
        self.canalyzing[i] = 1 - self.canalyzing[i]
        self.canalyzed[i] = 1 - self.canalyzed[i]


    def c_and_c_values_swapping(self, gene: int) -> None:
        """Swaps canalyzing and canalyzed values of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i1 = self.indices.index(gene)
        i2 = choice(list(set(range(self.arity)) - {i1}))  # picks another index different from <i1>
        self.canalyzing[i1], self.canalyzing[i2] = self.canalyzing[i2], self.canalyzing[i1]
        self.canalyzed[i1], self.canalyzed[i2] = self.canalyzed[i2], self.canalyzed[i1]
        self.indices[i1], self.indices[i2] = self.indices[i2], self.indices[i1]


    def c_and_c_values_insertion(self, gene: int) -> None:
        """Inserts update rule of given gene with randomly generated canalyzing and canalyzed values. <gene> can not
        occur in given NCF before.

        :param gene  id of gene"""

        i = randint(0, self.arity)
        self.canalyzing.insert(i, randint(0, 1))
        self.canalyzed.insert(i, randint(0, 1))
        self.indices.insert(i, gene)


    def c_and_c_values_removal(self, gene: int) -> None:
        """Removes canalyzing and canalyzed values of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.indices.index(gene)
        self.canalyzing.pop(i)
        self.canalyzed.pop(i)
        self.indices.pop(i)


def generate_rule(sign: Optional[bool]) -> Tuple[int, int]:
    """Generates random canalyzing and canalyzed pair depending on given sign of regulation

    :param sign  if True then chooses between AND, OR
                 if False then chooses between NAND, NOR
                 if None then chooses randomly of all 4
    :return tuple canalyzing, canalyzed value"""

    if sign is None:
        return choice([(1, 1), (0, 0), (0, 1), (1, 0)])
    return choice([(1, 1), (0, 0)] if sign else [(0, 1), (1, 0)])


class BooleanNetwork:
    """Class represents Boolean network (BN)

    Attributes:
        num_of_variables  number of network's Boolean variables
        functions         list of <num_of_variables> update functions in the form of NCFs"""

    def __init__(self, num_of_variables: int):
        self.num_of_variables: int = num_of_variables
        self.functions: List[UpdateFunction] = [UpdateFunction(i, num_of_variables) for i in range(num_of_variables)]

    @property
    def total_num_of_regulations(self):
        """Actual number of regulations in the entire network"""

        return sum(map(lambda x: x.arity, self.functions))


    def initialize_ncfs(self, regulators: Set[Regulation]) -> None:
        """Initializes NCF rules of given BN by using set of fixed regulations that can not be changed

        :param regulators  set of fixed regulations that can not be changed"""

        for regulation in regulators:
            self.add_regulator(regulation, True)


    def add_regulator(self, regulation: Regulation, fixed: bool, position: Optional[int] = None) -> None:
        """Adds new regulation to target gene specified in <regulation> on a specific position

        :param regulation  structured regulation object to be added
        :param fixed       True if regulator can not be removed afterwards, otherwise False
        :param position    specifies the position in NCF rule where to be inserted, if None, it is appended to
                           the last position in NCF rule, can not be greater than functions arity"""

        assert position is None or 0 <= position <= self.functions[regulation.target].arity
        self.functions[regulation.target].add_gene(regulation.source, fixed, position)
        canalyzing, canalyzed = generate_rule(regulation.sign)
        self.functions[regulation.target].add_rule(canalyzing, canalyzed, position)


    def mutate(self, num_of_genes: int, num_of_mutations: int) -> None:
        """Mutates particular number of genes of given BN and performs particular number of mutations of each gene

        :param num_of_genes      number of genes to be mutated
        :param num_of_mutations  number of mutations to be performed on each gene"""

        genes_to_mutate = choices(range(self.num_of_variables), k=num_of_genes)
        for gene in genes_to_mutate:
            self.functions[gene].mutate(self.num_of_variables, num_of_mutations, self.total_num_of_regulations)


    def to_aeon_string(self, perturbed_gene: int, isolated_variables: Set[int],
                       perturbation_state: Optional[bool] = None, ) -> str:
        """Returns string representation of given BN in .aeon format. Dimension of the model can be reduced due
        to the existence of isolated variables that will not be included in attractor analysis.

        :param perturbed_gene      denotes perturbed gene index, None if there is non perturbation
        :param isolated_variables  set of variables that do not regulate other variables and are not regulated too
        :param perturbation_state  False if perturbation is a knockout (KO) True if over-expression (OE),
                                   None if there is no perturbation
        :return AEON string representation of reduced model"""

        model_string = str()

        # adds all meaningful regulations to the model
        for i in range(self.num_of_variables):  # iterates over functions of given BN
            # perturbed gene acts as fixed input node therefore, its original regulations are skipped
            if i == perturbed_gene:
                continue

            for j in range(self.functions[i].arity):  # iterates over variables of i-th function
                reg_type = ">" if self.functions[i].canalyzing[j] == self.functions[i].canalyzed[j] else "|"
                model_string += "v_{} -{} v_{}\n".format(self.functions[i].indices[j], reg_type, i)

        # adds all meaningful logic update functions to the model
        for i in range(self.num_of_variables):

            if i not in isolated_variables:  # isolated variables are not added at all

                if i == perturbed_gene:  # perturbed gene has fixed update function - true (KO), false (OE)

                    if self.functions[i].indices != [i]:
                        # edge case: fixing variable that has only reflexive regulation results in its isolation
                        model_string += "$v_{}:{}\n".format(i, str(perturbation_state).lower())

                else:  # create update function iff gene is not perturbed and is not isolated
                    model_string += "{}\n".format(self.functions[i].to_aeon_string())

        return model_string


    def get_regulated_by(self, gene: int) -> Set[int]:
        """Returns set of genes that are regulated by given gene.

        :param gene  target regulator
        :return set of genes regulated by given gene"""

        result = set()
        for i in range(self.num_of_variables):
            for j in self.functions[i].indices:
                if gene == j:
                    result.add(i)

        return result


    def get_isolated_variables(self, perturbed_gene: int):
        """Returns set of isolated variables of given BN i.e., variables that have no regulators (inputs)
         and do not regulate any other variable (outputs) at the same time.

        :return set of isolated variables"""

        # perturbed gene is definitely input because it looses its update function so also its regulators
        input_variables = set() if perturbed_gene == -1 else {perturbed_gene}
        output_variables = set()

        # input variables
        for act_gene in range(self.num_of_variables):
            if not self.functions[act_gene].indices:  # only genes that have empty update function are inputs
                input_variables.add(act_gene)

        # output variables
        for act_gene in range(self.num_of_variables):
            regulates = self.get_regulated_by(act_gene)

            # gene is an output variable iff it regulates either no other genes or only perturbed gene (this regulation
            # will be formally removed in mutant type experiment by fixing value of perturbed gene)
            # in wild type experiment, <perturbed_gene> is equal to -1, which is an invalid gene index, so in wild
            # type experiment only genes that regulates no other genes are considered as output variables
            if regulates.issubset({perturbed_gene}):
                output_variables.add(act_gene)

        return output_variables.intersection(input_variables)


    def remove_perturbed_gene(self, perturbed_gene: int) -> bool:
        if perturbed_gene == -1:
            return False
        return len(self.get_regulated_by(perturbed_gene)) > 0


class Generation:
    """Class represents one generation of BNs of genetic algorithm

    Attributes:
        num_of_nets       number of networks in generation
        num_of_variables  number of gene variables in each network
        target_sinks      steady-states observed from input data, which are targeted by the inference
        networks          list of <num_of_nets> instances of <BooleanNetwork> that Generation consists of
        scores            fitness score for each network from <networks> in [0; 1]
                          describes how well the network fits the input data"""

    def __init__(self, num_of_nets: int, num_of_variables: int, target_sinks: Dict[int, List[Tuple[bool]]],
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
        :return instance of new generation created"""

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

            for j in sorted(self.target_sinks.keys()):  # iterates over perturbed gene indices of all experiments
                isolated_variables = self.networks[i].get_isolated_variables(j)
                # if WT then no gene is perturbed, if MT then looks at j-th variable of first (but basically any)
                # steady-state and derives type of perturbation (if 0 then it is KO, if 1 then it is OE)
                perturbation_type = None if j == -1 else self.target_sinks[j][0][j]
                aeon_model_string = self.networks[i].to_aeon_string(j, isolated_variables, perturbation_type)
                model, sag = get_symbolic_async_graph(aeon_model_string)
                observed_sinks = detect_steady_states(sag)
                total_score += evaluate_fitness(model, sag,
                                                reduce_attractors_dimension(self.target_sinks[j], isolated_variables),
                                                observed_sinks)
            self.scores[i] = total_score / len(self.target_sinks)  # final net's score is average score of all experiments


def create_initial_generation(num_of_nets: int, num_of_variables: int, input_constraints: Set[Regulation],
                              derived_constraints: Set[Regulation], sinks: Dict[int, List[State]]) -> Generation:
    """Creates the initial generation using regulations observed from input and derived constraints. Each BN is then
    mutated only once on one gene.

    :param num_of_nets          number of nets in generation
    :param num_of_variables     number of gene variables of each BN
    :param input_constraints    constraints specified in input file
    :param derived_constraints  constraints derived by gene-to-gene correlation
    :param sinks                target steady-state attractors specified in the input file, key corresponds to
                                the mutated gene, value contains the list of corresponding steady-state attractors
    :return instance of first generation"""

    init_gen = Generation(num_of_nets, num_of_variables, sinks)
    for net in init_gen.networks:
        net.initialize_ncfs(input_constraints.union(derived_constraints))
    init_gen.mutate(1, 1)  # allows only one mutation on one gene in the initial generation
    return init_gen


class BipartiteGraph:
    """Class represents instance of bipartite graph used for m-to-n comparison of BN attractors

    Attributes:
        target     vertices representing target steady-states (from input data)
        observed   vertices representing observed steady-states (from attractor analysis)
        distances  matrix m*n representing edges between each pair of target-observed steady-state
                   value distances[m][n] represent distance between m-th observed and n-th target steady-state"""

    def __init__(self, target_sinks, observed_sinks):
        self.target = target_sinks
        self.observed = observed_sinks
        # 2D array target x observed - distances[y][x] denotes manhattan dst between x-th target and y-th observed sink
        self.distances: List[List[Optional[int]]] = [[None] * len(self.target) for _ in range(len(self.observed))]
        self.calculate_distances()


    def calculate_distances(self):
        """Calculate distances between all pairs observed-target steady-state. Sets .distances matrix
        Distance between two states is equal to their Manhattan distance."""

        for i in range(len(self.observed)):
            for j in range(len(self.target)):
                self.distances[i][j] = manhattan_distance(self.observed[i], self.target[j])


    def minimal_weighted_assignment(self) -> Tuple[Optional[int], List[Tuple[int, int]]]:
        """Calculates minimal weighted assignment of given (possibly unbalanced) bipartite graph

        :return tuple - cost of minimal assignment, list of tuples of matching target-observed sink index pairs"""

        if not self.observed or not self.target:
            return None, []

        cost = np.array(self.distances)
        row_ind, col_ind = linear_sum_assignment(cost)
        return cost[row_ind, col_ind].sum(), list(zip(row_ind.tolist(), col_ind.tolist()))


def manhattan_distance(state1: Tuple[bool], state2: Tuple[bool]) -> int:
    """Calculates manhattan distance between two steady-states.

    :param state1  first attractor
    :param state2  second attractor
    :return their Manhattan distance"""

    assert len(state1) == len(state2), "State1:" + str(state1) + "State2: " + str(state2)
    result = 0
    for i in range(len(state1)):
        if state1[i] != state2[i]:
            result += 1
    return result


def evaluate_fitness(model, sag: SymbolicAsyncGraph, target_sinks: List[State], observed_sinks: List[State]) -> float:
    """Evaluates fitness of given BN depending on its steady-state attractors comparing to steady-states from
    given attractor data.

    :param model           biodivine_aeon.BooleanNetwork model of actual BN
    :param sag             Symbolic Asynchronous Graph of actual model
    :param target_sinks    steady-states of particular experiment from data
    :param observed_sinks  steady-states observed from model attractor analysis of particular experiment
    :return real number from [0;1] that determines BN's fitness"""

    bpg = BipartiteGraph(target_sinks, observed_sinks)
    cost, pairs = bpg.minimal_weighted_assignment()
    dimension = len(target_sinks[0])
    # weight of one variable in one state is computed as:
    # number of overlapping state tuples + number of overhanging states, multiplied by dimension of reduced model
    # therefore, each assigned tuple of states "behaves as one" and has weight equal to their reduced dimension and
    # each overhanging state alone has weight equal to its reduced dimension, one variable thus, have weight equal to
    # 1 / total number of variables where matching tuple of variables act as one
    matching_variables = 0
    total_variables = (abs(len(target_sinks) - len(observed_sinks)) + min(len(target_sinks), len(observed_sinks))) * dimension

    # if some states were matched, then total number of matching variables is equal to total number of variables
    # minus cost (variables that do not match), else it stays equal to 0 (initial value)
    if cost is not None:
        matching_variables += ((min(len(target_sinks), len(observed_sinks)) * dimension) - cost)

    # try to observe how many sinks absent and on which side
    if len(target_sinks) == len(observed_sinks):
        return matching_variables * total_variables

    elif len(target_sinks) > len(observed_sinks):
        # not ideal - some steady-states from data were not reached by given model,
        # check if missing steady-states are on some other type of attractor
        unmatched_states = get_unmatched_states(bpg, pairs, 0, len(target_sinks))
        for state in unmatched_states:
            if is_attractor_state(model, sag, state):
                matching_variables += dimension * 1/2  # penalty 1/3 for not being in single state

    else:  # len(target_sinks) < len(observed_sinks)
        # there is possibility that some steady-states were not caught while measuring, not a big problem if only few
        unmatched_states = get_unmatched_states(bpg, pairs, 1, len(observed_sinks))
        matching_variables += len(unmatched_states) * dimension * 3/4  # penalty for not being in data

    return matching_variables / total_variables


def get_unmatched_states(bpg: BipartiteGraph, matched_pairs: List[Tuple[int, int]],
                         position: int, total_number: int) -> List[State]:
    """Returns indices of states that were not matched in minimal weighted assignment.

    :param bpg            bipartite graph of steady-states of current BN
    :param matched_pairs  list of tuples of matched steady-states
    :param position       0 if some target steady-stated were not matched, 1 if observed were not matched
    :param total_number   total number of steady-states, specified by <position>
    :return list of states that were not matched in minimal weighted assignment"""

    assert len(bpg.target) != len(bpg.observed)
    matched = list(map(itemgetter(position), matched_pairs))
    unmatched_indices = set(range(total_number)) - set(matched)  # set difference returns unmatched indices
    unmatched_states = []
    overhung = bpg.target if position == 0 else bpg.observed  # link to either target of observed list of steady-states

    for i in unmatched_indices:
        unmatched_states.append(overhung[i])
    return unmatched_states


def reduce_attractors_dimension(sinks: List[State], isolated_variables: Set[int]) -> List[Tuple[Any]]:
    """Reduces the dimension of all entry states by eliminating states of isolated variables.

    :param sinks               entry list of steady-states
    :param isolated_variables  indices of variables that are isolated in actual model
    :return list of steady-states with reduced dimension"""

    return list(map(lambda s: tuple([s[i] for i in range(len(s)) if i not in isolated_variables]), sinks))


"""b = BooleanNetwork(4)
b.functions[0].indices.append(0)
b.functions[1].indices.extend([1])
b.functions[2].indices.extend([2])
b.functions[3].indices.append(3)
for i in range(-1, 4, 1):
    print("Experiment ", i)
    print(b.get_isolated_variables(i))
"""