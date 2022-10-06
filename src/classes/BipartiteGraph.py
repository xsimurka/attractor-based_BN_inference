from typing import Tuple, Optional, List
import numpy as np
from scipy.optimize import linear_sum_assignment
import src.utils as utils


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
                self.distances[i][j] = utils.manhattan_distance(self.observed[i], self.target[j])

    def minimal_weighted_assignment(self) -> Tuple[Optional[int], List[Tuple[int, int]]]:
        """Calculates minimal weighted assignment of given (possibly unbalanced) bipartite graph

        :return tuple  cost of minimal assignment, list of tuples of matching target-observed sink index pairs"""

        if not self.observed or not self.target:
            return None, []

        cost = np.array(self.distances)
        row_ind, col_ind = linear_sum_assignment(cost)
        return cost[row_ind, col_ind].sum(), list(zip(row_ind.tolist(), col_ind.tolist()))
