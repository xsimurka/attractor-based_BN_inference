from typing import Tuple, Dict, List, Set, Optional
import src.classes.Regulation as reg
import pandas as pd
from operator import itemgetter
import itertools
from random import random

State = Tuple[bool]


class BNInfo:
    """ Class holds information about the target network.

    Attributes:
        num_of_vars       number of variables of target network
        input_genes       set of input nodes i.e., nodes that are not regulated by other genes
        output_genes      set of output nodes i.e., nodes that do not regulate other genes
        wt_sinks          sinks that corresponds to wild-type experiment
        ko_sinks          dictionary where key corresponds to perturbed gene and value to list of steady-states of
                          its knock-out experiment
        oe_sinks          dictionary where key corresponds to perturbed gene and value to list of steady-states of
                          its over-expression experiment"""

    def __init__(self, num_of_vars: int, input_genes: Set[int], output_genes: Set[int]):
        self.num_of_vars = num_of_vars
        self.input_genes = input_genes
        self.output_genes = output_genes
        self.wt_sinks: List[State] = []
        self.ko_sinks: Dict[int, List[State]] = {}
        self.oe_sinks: Dict[int, List[State]] = {}

    def derive_constraints(self, threshold: float) -> Set[reg.Regulation]:
        """Computes Pearson's correlation between each pair of genes, if the absolute value is greater than given
        <threshold> then suggests such regulation to be added later during the inference.

        :param threshold  the value that the correlation must exceed in order to be suggested
        :return           list of suggested undirected regulations"""

        print("Deriving constraints from gene-to-gene correlations...", end="")
        derived_regulations: Set[reg.Regulation] = set()
        for i in range(self.num_of_vars):
            for j in range(i + 1, self.num_of_vars, 1):
                corr = self.calculate_correlation(i, j)
                if abs(corr) >= threshold:
                    r = self.create_regulation(i, j, corr)
                    if r is not None:
                        derived_regulations.add(r)
        print(" done.")

        if derived_regulations:
            print("Derived non-direct regulations:")
            for regulation in derived_regulations:
                print("v_{} {}={} v_{}".format(regulation.source,
                                               "<" if regulation.sign else "|",
                                               ">" if regulation.sign else "|",
                                               regulation.target))
        else:
            print("No derived regulations.")

        return derived_regulations

    def create_regulation(self, gene1: int, gene2: int, correlation: float) -> Optional[reg.Regulation]:
        """

        """

        if gene1 in self.input_genes and gene2 in self.input_genes:
            return None
        if gene1 in self.output_genes and gene2 in self.output_genes:
            return None
        if gene1 in self.input_genes and gene2 in self.output_genes:
            return reg.Regulation(gene1, gene2, correlation > 0)
        if gene2 in self.input_genes and gene1 in self.output_genes:
            return reg.Regulation(gene2, gene1, correlation > 0)
        if random() < 0.5:
            gene1, gene2 = gene2, gene1
        return reg.Regulation(gene1, gene2, correlation > 0)

    def calculate_correlation(self, gene1: int, gene2: int) -> float:
        """Calculates single correlation between given tuple of genes.

        :param gene1  first gene index
        :param gene2  second gene index
        :return       correlation between given tuple of genes"""

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

        ko_steady_states = list(itertools.chain.from_iterable(self.ko_sinks.values()))
        return list(map(itemgetter(gene), ko_steady_states))

    def get_gene_oe_values(self, gene: int) -> List[bool]:
        """Returns vector of <gene> over-expression values.

        :param gene  index of gene
        :return vector of its over-expression boolean values"""

        oe_steady_states = list(itertools.chain.from_iterable(self.oe_sinks.values()))
        return list(map(itemgetter(gene), oe_steady_states))
