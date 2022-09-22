from biodivine_aeon import *
from typing import List, Tuple


State = Tuple[bool]


def get_symbolic_async_graph(model_string) -> Tuple[BooleanNetwork, SymbolicAsyncGraph]:
	"""Precomputes model computational structures from string representation, which is essential for iterative approach
	:param model_string  string representation of Boolean model
	:return tuple of computational structures of given model"""

	model = BooleanNetwork.from_aeon(model_string)
	return model, SymbolicAsyncGraph(model)


def detect_steady_states(sag: SymbolicAsyncGraph) -> List[State]:
	"""Returns list of steady states of given model represented by Symbolic Asynchronous graph.
	:param sag  Symbolic Asynchronous graph of Boolean model
	:return list of steady states"""

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
