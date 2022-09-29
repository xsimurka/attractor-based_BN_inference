from math import sqrt
from typing import List, Tuple, Dict, Set
from src.classes.Regulation import Regulation
import src.classes.BNInfo as bn
import pandas as pd
from operator import itemgetter

State = Tuple[bool]


def derive_constraints(sinks: bn.BNInfo, threshold: float) -> Set[Regulation]:
    derived_regulations: Set[Regulation] = set()
    for i in range(sinks.num_of_variables):
        for j in range(i, sinks.num_of_variables, 1):
            corr = calculate_correlation(sinks, i, j)
            if abs(corr) < threshold:
                continue
            coh_ij = calculate_gene_to_gene_coherency(i, j)
            coh_ji = calculate_gene_to_gene_coherency(j, i)
            if corr * coh_ij > 0:
                derived_regulations.add(Regulation(i, j, corr > 0))
            if corr * coh_ji > 0:
                derived_regulations.add(Regulation(j, i, corr > 0))
    return derived_regulations


def calculate_correlation(sinks: bn.BNInfo, gene1: int, gene2: int) -> float:
    gene1_values = pd.DataFrame(get_gene_sinks_values(sinks, gene1))
    gene2_values = pd.DataFrame(get_gene_sinks_values(sinks, gene2))
    return float(gene1_values.corrwith(gene2_values)[0])


def calculate_gene_to_gene_coherency(gene1: int, gene2: int) -> float:
    pass


def get_gene_sinks_values(bn_info: bn.BNInfo, gene: int) -> List[bool]:
    gene_values = []
    gene_values.extend(list(map(itemgetter(gene), bn_info.wt_sinks)))
    gene_values.extend(list(map(itemgetter(gene), bn_info.ko_sinks.values())))
    gene_values.extend(list(map(itemgetter(gene), bn_info.oe_sinks.values())))
    return gene_values


print(pd.DataFrame(map(bool, [1, 1, 1, 0, 1, 0, 1, 1])).corrwith(pd.DataFrame(map(bool, [0, 1, 0, 0, 0, 1, 1, 0])), 0,
                                                                 False, "pearson"))
print(float(
    pd.DataFrame([0, 1, 0, 0, 0, 1, 1, 0]).corrwith(pd.DataFrame([1, 1, 1, 0, 1, 0, 1, 1]), 0, False, "pearson")[0]))
