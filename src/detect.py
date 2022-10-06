# Author: RNDr. Samuel Pastva, Ph.D.


from biodivine_aeon import *
from typing import List, Tuple

State = Tuple[bool]


def get_model_computational_structures(model_string) -> Tuple[BooleanNetwork, SymbolicAsyncGraph]:
    """Precomputes model computational structures from string representation, which is essential for iterative approach
    :param model_string  string representation of Boolean model
    :return              tuple of computational structures of given model"""

    model = BooleanNetwork.from_aeon(model_string)
    return model, SymbolicAsyncGraph(model)


def detect_steady_states(sag: SymbolicAsyncGraph) -> List[State]:
    """Returns list of steady states of given model represented by Symbolic Asynchronous graph.
    :param sag  Symbolic Asynchronous graph of Boolean model
    :return     list of steady states"""

    candidates = []
    for var in sag.network().variables():
        can_post = sag.var_can_post(var, sag.unit_colored_vertices())
        candidates.append(sag.unit_colored_vertices().minus(can_post))

    while len(candidates) > 1:
        candidates.sort(reverse=True, key=lambda x: x.symbolic_size())
        item = candidates.pop()
        candidates = [item.intersect(x) for x in candidates]

    sinks = candidates[0].vertices().list_vertices()
    if not sinks:
        return []
    return list(map(tuple, sinks))


def is_attractor_state(model, sag: SymbolicAsyncGraph, state: Tuple[bool]) -> bool:
    """Finds out whether given states lies in any kind of attractor of given model. For the iterative approach
    it is essential to compute model and symbolic graph just once as it is more time demanding.

    :param model  biodivine_aeon.BooleanNetwork model of actual network
    :param sag    Symbolic Asynchronous Graph of actual network
    :param state  desired state to check
    :return True if state lies in attractor of given model, otherwise False"""

    init_state = sag.fix_vertex(list(state))
    bwd = init_state
    while True:
        done = True
        for var in reversed(model.variables()):
            step = sag.var_pre(var, bwd)
            if not step.is_subset(bwd):
                bwd = bwd.union(step)
                done = False
                break
        if done:
            break

    attractor = init_state
    while True:
        done = True
        for var in reversed(model.variables()):
            step = sag.var_post(var, attractor)

            if not step.is_subset(attractor):
                bad_colors = step.minus(bwd).colors()
                attractor = attractor.union(step).minus_colors(bad_colors)
                done = False
                break
        if done:
            break

    return len(attractor.vertices().list_vertices()) > 0
