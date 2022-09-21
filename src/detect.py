from biodivine_aeon import *
from typing import List, Tuple


State = Tuple[bool]


def get_symbolic_async_graph(model_string) -> Tuple[BooleanNetwork, SymbolicAsyncGraph]:
	"""Docstring TODO"""
	model = BooleanNetwork.from_aeon(model_string)
	return model, SymbolicAsyncGraph(model)


def detect_steady_states(sag: SymbolicAsyncGraph) -> List[State]:
	"""Docstring TODO"""

	candidates = []
	for var in sag.network().variables():
		can_post = sag.var_can_post(var, sag.unit_colored_vertices())
		candidates.append(sag.unit_colored_vertices().minus(can_post))

	while len(candidates) > 1:
		candidates.sort(reverse=True, key=lambda x: x.symbolic_size())
		item = candidates.pop()
		candidates = [item.intersect(x) for x in candidates]

	if not candidates[0].vertices().list_vertices():
		return []
	return list(map(tuple, candidates[0].vertices().list_vertices()))
