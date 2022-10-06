from typing import Tuple, Dict, List
import src.classes.Regulation as reg
import pandas as pd
from operator import itemgetter

State = Tuple[bool]


class BNInfo:
    """ Class holds information about the input steady-state data
    Attributes:
        num_of_variables  number of variables of target network
        wt_sinks          sinks that corresponds to wild-type experiment
        ko_sinks          dictionary where key corresponds to perturbed gene and value to list of steady-states of
                          its knock-out experiment
        oe_sinks          dictionary where key corresponds to perturbed gene and value to list of steady-states of
                          its over-expression experiment"""

    def __init__(self, num_of_variables: int):
        self.num_of_variables = num_of_variables
        self.wt_sinks: List[State] = []
        self.ko_sinks: Dict[int, List[State]] = {}
        self.oe_sinks: Dict[int, List[State]] = {}

    def derive_constraints(self, threshold: float) -> List[reg.Regulation]:
        """Computes Pearson's correlation between each pair of genes, if the absolute value is greater than given
        <threshold> then suggests such regulation to be added later during the inference.

        :param threshold  the value that the correlation must exceed in order to be suggested
        :return          list of suggested undirected regulations"""

        print("Deriving constraints from gene-to-gene correlations...")
        derived_regulations: List[reg.Regulation] = list()
        for i in range(self.num_of_variables):
            for j in range(i + 1, self.num_of_variables, 1):
                corr = self.calculate_correlation(i, j)
                if abs(corr) >= threshold:
                    derived_regulations.append(reg.Regulation(i, j, corr > 0, False))
        print(" done.")

        if derived_regulations:
            print("Derived non-direct regulations:")
            for regulation in derived_regulations:
                print("v_{} - v_{} sign: {}".format(regulation.source, regulation.target, regulation.sign))
        else:
            print("No derived regulations.")

        return derived_regulations

    def calculate_correlation(self, gene1: int, gene2: int) -> float:
        """Calculates single correlation between given tuple of genes.

        :param gene1  first gene index
        :param gene2  second gene index
        :return       correlation between given tuple of genes """

        gene1_values = pd.DataFrame(self.get_gene_sinks_values(gene1))
        gene2_values = pd.DataFrame(self.get_gene_sinks_values(gene2))
        return float(gene1_values.corrwith(gene2_values)[0])

    def get_gene_sinks_values(self, gene: int) -> List[bool]:
        """Returns vector of boolean values of given <gene> over all the input experiments.

        :param gene  index of gene
        :return vector of its boolean values"""

        gene_values = []
        gene_values.extend(self.get_gene_wt_values(gene))
        gene_values.extend(self.get_gene_ko_values(gene))
        gene_values.extend(self.get_gene_oe_values(gene))
        return gene_values

    def get_gene_wt_values(self, gene: int) -> List[bool]:
        """Returns vector of <gene> wild-type values.

        :param gene  index of gene
        :return vector of its wild-type boolean values"""

        return list(map(itemgetter(gene), self.wt_sinks))

    def get_gene_ko_values(self, gene: int) -> List[bool]:
        """Returns vector of <gene> knock-out values.

        :param gene  index of gene
        :return vector of its knock-out boolean values"""

        return list(map(itemgetter(gene), self.ko_sinks.values()))

    def get_gene_oe_values(self, gene: int) -> List[bool]:
        """Returns vector of <gene> over-expression values.

        :param gene  index of gene
        :return vector of its over-expression boolean values"""

        return list(map(itemgetter(gene), self.oe_sinks.values()))
