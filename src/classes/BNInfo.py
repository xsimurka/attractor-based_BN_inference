from typing import Tuple, Dict, List
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

    def derive_constraints(self, threshold: float) -> List[reg.Regulation]:
        derived_regulations: List[reg.Regulation] = list()
        for i in range(self.num_of_variables):
            for j in range(i + 1, self.num_of_variables, 1):
                corr = self.calculate_correlation(i, j)
                if abs(corr) >= threshold:
                    derived_regulations.append(reg.Regulation(i, j, corr > 0, False))
        return derived_regulations

    def calculate_correlation(self, gene1: int, gene2: int) -> float:
        gene1_values = pd.DataFrame(self.get_gene_sinks_values(gene1))
        gene2_values = pd.DataFrame(self.get_gene_sinks_values(gene2))
        return float(gene1_values.corrwith(gene2_values)[0])

    def get_gene_sinks_values(self, gene: int) -> List[bool]:
        gene_values = []
        gene_values.extend(self.get_gene_wt_values(gene))
        gene_values.extend(self.get_gene_ko_values(gene))
        gene_values.extend(self.get_gene_oe_values(gene))
        return gene_values

    def get_gene_wt_values(self, gene: int) -> List[bool]:
        return list(map(itemgetter(gene), self.wt_sinks))

    def get_gene_ko_values(self, gene: int) -> List[bool]:
        return list(map(itemgetter(gene), self.ko_sinks.values()))

    def get_gene_oe_values(self, gene: int) -> List[bool]:
        return list(map(itemgetter(gene), self.oe_sinks.values()))
