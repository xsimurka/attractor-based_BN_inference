from copy import deepcopy
from random import randint, choices, choice
from typing import Tuple, Optional, List
import src.classes.BNInfo as bn


class UpdateFunction:
    """Class represents boolean update rule in the form of Nested Canalyzing Function (NCF)

    Attributes:
        gene            id of regulated gene
        target_bn_info  object that holds the information about target BN
        canalyzing      canalyzing values
        canalyzed       canalyzed values
        regulators      ids of genes that regulates <gene>
        fixed           set of fixed regulators of actual gene"""

    def __init__(self, gene: int, target_bn_info: bn.BNInfo):
        self.gene = gene
        self.target_bn_info = target_bn_info
        self.canalyzing = []
        self.canalyzed = []
        self.regulators = []
        self.fixed = set()

    def deepcopy__(self):
        result: UpdateFunction = UpdateFunction(self.gene, self.target_bn_info)
        result.canalyzed = deepcopy(self.canalyzed)
        result.canalyzing = deepcopy(self.canalyzing)
        result.regulators = deepcopy(self.regulators)
        result.fixed = deepcopy(self.fixed)
        return result

    @property
    def non_fixed(self) -> int:
        """Number of non fixed regulators"""
        return len(self.regulators) - len(self.fixed)

    @property
    def arity(self) -> int:
        """Actual number of regulators"""
        return len(self.regulators)

    def to_aeon_string(self) -> str:
        """Returns string representation of given NCF in .aeon format

        :return  aeon string representation"""

        if self.arity == 0:  # an empty NCF formula can not occur in AEON string (results in an error)
            return str()

        # the most nested regulator will not have brackets around
        function_string = "{}v_{}".format("" if self.canalyzing[-1] == self.canalyzed[-1] else "!", self.regulators[-1])

        for i in reversed(range(self.arity - 1)):  # without the last regulator already
            logic_operator = "&" if self.canalyzed[i] == 0 else "|"
            parity = "" if self.canalyzing[i] == self.canalyzed[i] else "!"
            function_string = "({}v_{} {} ".format(parity, self.regulators[i], logic_operator) + function_string + ")"

        return "$v_{}:{}".format(self.gene, function_string)

    def add_gene(self, gene: int, fixed: bool, position: Optional[int] = None) -> None:
        """Adds new gene index to list of regulators of given NCF on given position. For time complexity reason it is
        differentiated between appending on the last position and inserting elsewhere.

        :param gene      id of added gene
        :param fixed     True if regulator can not be removed later during the inference, otherwise False
        :param position  specifies the position in NCF rule where to be inserted
                         if None, it is appended to the last position in NCF rule"""

        if position is None:
            self.regulators.append(gene)
        else:
            self.regulators.insert(position, gene)
        if fixed:
            self.fixed.add(gene)

    def add_rule(self, canalyzing: int, canalyzed: int, position: Optional[int] = None) -> None:
        """Adds given canalyzing and canalyzed value to given NCF rule on given position. For time complexity reason
        it is differentiated between appending on the last position and inserting elsewhere.

        :param canalyzing  canalyzing value
        :param canalyzed   canalyzed value
        :param position    specifies the position in NCF rule where to be inserted
                           if None, it is appended to the last position in NCF rule"""

        if position is None:
            self.canalyzing.append(canalyzing)
            self.canalyzed.append(canalyzed)
        else:
            self.canalyzing.insert(position, canalyzing)
            self.canalyzed.insert(position, canalyzed)

    def evaluate(self, interpretation: Tuple[int]) -> bool:
        """Returns the valuation of given NCF rule of given interpretation
        If variable from interpretation matches its corresponding canalyzing value, function returns corresponding
        canalyzed value. If no variable matches canalyzing value, functions evaluates itself to negation of the last
        canalyzed value.

        :param interpretation  Boolean vector specifying interpretation of atomic variables
        :return                Boolean valuation of given NCF using given interpretation"""

        assert (len(interpretation) == len(self.canalyzing) == len(self.canalyzed))
        for i in range(len(interpretation)):
            if interpretation[i] == self.canalyzing[i]:
                return bool(self.canalyzed[i])

        return not bool(self.canalyzed[-1])

    def select_reasonable_mutations(self, gene: int, total_num_of_regulations: int) -> List[str]:
        """Selects all reasonable mutations for particular selected gene according to its occurrence in NCF rule.
        Leaves possibility that no mutation is available - if all regulations of picked gene are fixed and model is
        already too dense.

        :param gene                      id of regulator to be checked
        :param total_num_of_regulations  the number of all current regulations in the network
        :return                          list of possible mutations that are reasonable to be done on given gene"""

        mutations = []
        if gene in self.regulators:
            if self.arity >= 2:
                mutations.append("c_and_c_values_swapping")

            if gene not in self.fixed:
                mutations.append("canalyzing_value_reversion")
                mutations.append("canalyzed_value_reversion")
                mutations.append("c_and_c_values_removal")

            if self.get_c_and_c_values(gene) in {(0, 1), (1, 0)}:
                mutations.append("c_and_c_values_reversion")

        elif total_num_of_regulations < (
                self.target_bn_info.num_of_vars ** 2) / 4:  # restricted number of regulations in network
            mutations.append("c_and_c_values_insertion")

        return mutations

    def mutate(self, num_of_mutations: int, total_num_of_regulations: int) -> None:
        """Selects randomly given amount of regulators to be mutated, select possible mutations of each
        regulators, performs chosen mutations

        :param num_of_mutations          number of mutations to be done
        :param total_num_of_regulations  total number of regulations in the whole network"""

        non_output_genes = list(set(range(self.target_bn_info.num_of_vars)) - self.target_bn_info.output_genes)
        if 5 in non_output_genes:
            pass
        regulators_to_mutate = choices(non_output_genes, k=num_of_mutations)

        for regulator in regulators_to_mutate:
            mutations = self.select_reasonable_mutations(regulator, total_num_of_regulations)
            if mutations:
                getattr(self, choice(mutations))(regulator)

    def get_c_and_c_values(self, gene: int) -> Optional[Tuple[int, int]]:
        """Returns canalyzing and canalyzed value of particular gene id if occurred in given NCF rule

        :param gene  id of gene
        :return      tuple with canalyzing and canalyzed value if gene occurs in given NCF rule, otherwise None"""

        try:
            i = self.regulators.index(gene)
        except ValueError:
            return None
        return self.canalyzing[i], self.canalyzed[i]

    def get_nth_non_fixed_gene_index(self, n: int) -> int:
        """Returns id of nth non fixed regulator of given NCF rule

        :param n
        :return  index of nth non fixed regulator, -1 if less than n regulators are non fixed"""

        result = -1
        if n > self.non_fixed:
            return result

        while n > 0:
            result += 1
            if self.regulators[result] not in self.fixed:
                n -= 1

        return result

    def canalyzing_value_reversion(self, gene: int) -> None:
        """Reverts canalyzing value of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.regulators.index(gene)
        self.canalyzing[i] = 1 - self.canalyzing[i]

    def canalyzed_value_reversion(self, gene: int) -> None:
        """Reverts canalyzed value of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.regulators.index(gene)
        self.canalyzed[i] = 1 - self.canalyzed[i]

    def c_and_c_values_reversion(self, gene: int) -> None:
        """Reverts canalyzing and canalyzed values of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.regulators.index(gene)
        self.canalyzing[i] = 1 - self.canalyzing[i]
        self.canalyzed[i] = 1 - self.canalyzed[i]

    def c_and_c_values_swapping(self, gene: int) -> None:
        """Swaps canalyzing and canalyzed values of given gene. Given NCF has to contain <gene> and one more regulator.

        :param gene  id of gene"""

        i1 = self.regulators.index(gene)
        i2 = choice(list(set(range(self.arity)) - {i1}))  # picks another index different from <i1>
        self.canalyzing[i1], self.canalyzing[i2] = self.canalyzing[i2], self.canalyzing[i1]
        self.canalyzed[i1], self.canalyzed[i2] = self.canalyzed[i2], self.canalyzed[i1]
        self.regulators[i1], self.regulators[i2] = self.regulators[i2], self.regulators[i1]

    def c_and_c_values_insertion(self, gene: int) -> None:
        """Inserts update rule of given gene with randomly generated canalyzing and canalyzed values. <gene> can not
        occur in given NCF before.

        :param gene  id of gene"""

        i = randint(0, self.arity)
        self.canalyzing.insert(i, randint(0, 1))
        self.canalyzed.insert(i, randint(0, 1))
        self.regulators.insert(i, gene)

    def c_and_c_values_removal(self, gene: int) -> None:
        """Removes canalyzing and canalyzed values of given gene. Given NCF has to contain <gene>.

        :param gene  id of gene"""

        i = self.regulators.index(gene)
        self.canalyzing.pop(i)
        self.canalyzed.pop(i)
        self.regulators.pop(i)
