from biodivine_aeon import *
from typing import Tuple


def is_attractor_state(model, sag: SymbolicAsyncGraph, state: Tuple[bool]):
	"""Docstring TODO"""

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
