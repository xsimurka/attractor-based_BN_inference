from biodivine_aeon import *
from typing import Tuple


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
