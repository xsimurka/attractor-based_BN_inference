from random import choice
from biodivine_aeon import SymbolicAsyncGraph
from typing import Set, Tuple, Optional, List, Any
from operator import itemgetter
from detect import is_attractor_state
from parse_input import BNInfo
import src.classes.Regulation as reg
import src.classes.Generation as gen
import src.classes.BipartiteGraph as bg


State = Tuple[bool]


def generate_rule(sign: Optional[bool]) -> Tuple[int, int]:
    """Generates random canalyzing and canalyzed pair depending on given sign of regulation

    :param sign  if True then chooses between AND, OR
                 if False then chooses between NAND, NOR
                 if None then chooses randomly of all 4
    :return      tuple canalyzing, canalyzed value"""

    if sign is None:
        return choice([(1, 1), (0, 0), (0, 1), (1, 0)])
    return choice([(1, 1), (0, 0)] if sign else [(0, 1), (1, 0)])



def create_initial_generation(num_of_nets: int, num_of_variables: int, input_constraints: Set[reg.Regulation],
                              derived_constraints: Set[reg.Regulation], sinks: BNInfo) -> gen.Generation:
    """Creates the initial generation using regulations observed from input and derived constraints. Each BN is then
    mutated only once on one gene.

    :param num_of_nets          number of nets in generation
    :param num_of_variables     number of gene variables of each BN
    :param input_constraints    constraints specified in input file
    :param derived_constraints  constraints derived by gene-to-gene correlation
    :param sinks                target steady-state attractors specified in the input file, key corresponds to
                                the mutated gene, value contains the list of corresponding steady-state attractors
    :return                     instance of first generation"""

    init_gen = gen.Generation(num_of_nets, num_of_variables, sinks)
    for net in init_gen.networks:
        net.initialize_ncfs(input_constraints.union(derived_constraints))
    init_gen.mutate(1, 1)  # allows only one mutation on one gene in the initial generation
    return init_gen



def manhattan_distance(state1: Tuple[bool], state2: Tuple[bool]) -> int:
    """Calculates manhattan distance between two steady-states.

    :param state1  first attractor
    :param state2  second attractor
    :return        their Manhattan distance"""

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
    :return                real number from [0;1] that determines BN's fitness"""

    bpg = bg.BipartiteGraph(target_sinks, observed_sinks)
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
    if len(target_sinks) > len(observed_sinks):
        # not ideal - some steady-states from data were not reached by given model,
        # check if missing steady-states are on some other type of attractor
        unmatched_states = get_unmatched_states(bpg, pairs, 0, len(target_sinks))
        for state in unmatched_states:
            if is_attractor_state(model, sag, state):
                matching_variables += dimension * 1/2  # penalty 1/3 for not being in single state

    elif len(target_sinks) < len(observed_sinks):
        # there is possibility that some steady-states were not caught while measuring, not a big problem if only few
        unmatched_states = get_unmatched_states(bpg, pairs, 1, len(observed_sinks))
        matching_variables += len(unmatched_states) * dimension * 3/4  # penalty for not being in data

    # if target == observed then no correction is needed
    return matching_variables / total_variables


def get_unmatched_states(bpg: bg.BipartiteGraph, matched_pairs: List[Tuple[int, int]],
                         position: int, total_number: int) -> List[State]:
    """Returns indices of states that were not matched in minimal weighted assignment.

    :param bpg            bipartite graph of steady-states of current BN
    :param matched_pairs  list of tuples of matched steady-states
    :param position       0 if some target steady-stated were not matched, 1 if observed were not matched
    :param total_number   total number of steady-states, specified by <position>
    :return               list of states that were not matched in minimal weighted assignment"""

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
    :return                    list of steady-states with reduced dimension"""

    return list(map(lambda s: tuple([s[i] for i in range(len(s)) if i not in isolated_variables]), sinks))
