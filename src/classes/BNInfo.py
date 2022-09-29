from typing import Tuple, Dict, List, Set
import src.classes.Regulation as reg
import pandas as pd
from operator import itemgetter

State = Tuple[bool]


class BNInfo:

    def __init__(self, num_of_variables: int):
        self.num_of_variables = num_of_variables
        self.wt_sinks: List[State] = []
        self.ko_sinks: Dict[int, List[State]] = {}
        self.oe_sinks: Dict[int, List[State]] = {}

    def derive_constraints(self, threshold: float) -> Set[reg.Regulation]:
        derived_regulations: Set[reg.Regulation] = set()
        for i in range(self.num_of_variables):
            for j in range(i, self.num_of_variables, 1):
                corr = self.calculate_correlation(i, j)
                if abs(corr) < threshold:
                    continue
                coh_ij = self.calculate_gene_to_gene_coherency(i, j)
                coh_ji = self.calculate_gene_to_gene_coherency(j, i)
                if corr * coh_ij > 0:
                    derived_regulations.add(reg.Regulation(i, j, corr > 0))
                if corr * coh_ji > 0:
                    derived_regulations.add(reg.Regulation(j, i, corr > 0))
        return derived_regulations

    def calculate_correlation(self, gene1: int, gene2: int) -> float:
        gene1_values = pd.DataFrame(self.get_gene_sinks_values(gene1))
        gene2_values = pd.DataFrame(self.get_gene_sinks_values(gene2))
        return float(gene1_values.corrwith(gene2_values)[0])

    def calculate_gene_to_gene_coherency(self, gene1: int, gene2: int) -> float:
        pass

    def get_gene_sinks_values(self, gene: int) -> List[bool]:
        gene_values = []
        gene_values.extend(list(map(itemgetter(gene), self.wt_sinks)))
        gene_values.extend(list(map(itemgetter(gene), self.ko_sinks.values())))
        gene_values.extend(list(map(itemgetter(gene), self.oe_sinks.values())))
        return gene_values
