from random import choice
from typing import Set, Tuple, Optional, List, Any
from operator import itemgetter
import src.classes.BNInfo as bn
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
                              derived_constraints: Set[reg.Regulation], sinks: bn.BNInfo) -> gen.Generation:
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
