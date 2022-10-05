from typing import Set, Optional, List
from random import choices
from src.classes.UpdateFunction import UpdateFunction
from src.classes.Regulation import Regulation
import src.utils as utils


class BooleanNetwork:
    """Class represents Boolean network (BN)

    Attributes:
        num_of_variables  number of network's Boolean variables
        functions         list of <num_of_variables> update functions in the form of NCFs"""

    def __init__(self, num_of_variables: int):
        self.num_of_variables: int = num_of_variables
        self.functions: List[UpdateFunction] = [UpdateFunction(i, num_of_variables) for i in range(num_of_variables)]

    @property
    def total_num_of_regulations(self):
        """Actual number of regulations in the entire network"""

        return sum(map(lambda x: x.arity, self.functions))

    def initialize_ncfs(self, regulators: List[Regulation]) -> None:
        """Initializes NCF rules of given BN by using set of fixed regulations that can not be changed

        :param regulators  set of fixed regulations that can not be changed"""

        for regulation in regulators:
            self.add_regulator(regulation, True)

    def add_regulator(self, regulation: Regulation, fixed: bool, position: Optional[int] = None) -> None:
        """Adds new regulation to target gene specified in <regulation> on a specific position

        :param regulation  structured regulation object to be added
        :param fixed       True if regulator can not be removed afterwards, otherwise False
        :param position    specifies the position in NCF rule where to be inserted, if None, it is appended to
                           the last position in NCF rule, can not be greater than functions arity"""

        assert position is None or 0 <= position <= self.functions[regulation.target].arity
        self.functions[regulation.target].add_gene(regulation.source, fixed, position)
        canalyzing, canalyzed = utils.generate_rule(regulation.sign)
        self.functions[regulation.target].add_rule(canalyzing, canalyzed, position)

    def mutate(self, num_of_genes: int, num_of_mutations: int) -> None:
        """Mutates particular number of genes of given BN and performs particular number of mutations of each gene

        :param num_of_genes      number of genes to be mutated
        :param num_of_mutations  number of mutations to be performed on each gene"""

        genes_to_mutate = choices(range(self.num_of_variables), k=num_of_genes)
        for gene in genes_to_mutate:
            self.functions[gene].mutate(self.num_of_variables, num_of_mutations, self.total_num_of_regulations)

    def to_aeon_string(self, perturbed_gene: int, isolated_variables: Set[int],
                       perturbation_state: Optional[bool] = None, ) -> str:
        """Returns string representation of given BN in .aeon format. Dimension of the model can be reduced due
        to the existence of isolated variables that will not be included in attractor analysis.

        :param perturbed_gene      denotes perturbed gene index, None if there is non perturbation
        :param isolated_variables  set of variables that do not regulate other variables and are not regulated too
        :param perturbation_state  False if perturbation is a knockout (KO) True if over-expression (OE),
                                   None if there is no perturbation
        :return                    AEON string representation of reduced model"""

        model_string = str()

        # adds all meaningful regulations to the model
        for i in range(self.num_of_variables):  # iterates over functions of given BN
            # perturbed gene acts as fixed input node therefore, its original regulations are skipped
            if i == perturbed_gene:
                continue

            for j in range(self.functions[i].arity):  # iterates over variables of i-th function
                reg_type = ">" if self.functions[i].canalyzing[j] == self.functions[i].canalyzed[j] else "|"
                model_string += "v_{} -{} v_{}\n".format(self.functions[i].indices[j], reg_type, i)

        # adds all meaningful logic update functions to the model
        for i in range(self.num_of_variables):

            if i not in isolated_variables:  # isolated variables are not added at all

                if i == perturbed_gene:  # perturbed gene has fixed update function - true (KO), false (OE)

                    if self.functions[i].indices != [i]:
                        # edge case: fixing variable that has only reflexive regulation results in its isolation
                        model_string += "$v_{}:{}\n".format(i, str(perturbation_state).lower())

                else:  # create update function iff gene is not perturbed and is not isolated
                    model_string += "{}\n".format(self.functions[i].to_aeon_string())

        return model_string

    def get_regulated_by(self, gene: int) -> Set[int]:
        """Returns set of genes that are regulated by given gene.

        :param gene  target regulator
        :return      set of genes regulated by given gene"""

        result = set()
        for i in range(self.num_of_variables):
            for j in self.functions[i].indices:
                if gene == j:
                    result.add(i)

        return result

    def get_isolated_variables(self, perturbed_gene: int):
        """Returns set of isolated variables of given BN i.e., variables that have no regulators (inputs)
         and do not regulate any other variable (outputs) at the same time.

        :return  set of isolated variables"""

        # perturbed gene is definitely input because it looses its update function so also its regulators
        input_variables = set() if perturbed_gene == -1 else {perturbed_gene}
        output_variables = set()

        # input variables
        for act_gene in range(self.num_of_variables):
            if not self.functions[act_gene].indices:  # only genes that have empty update function are inputs
                input_variables.add(act_gene)

        # output variables
        for act_gene in range(self.num_of_variables):
            regulates = self.get_regulated_by(act_gene)

            # gene is an output variable iff it regulates either no other genes or only perturbed gene (this regulation
            # will be formally removed in mutant type experiment by fixing value of perturbed gene)
            # in wild type experiment, <perturbed_gene> is equal to -1, which is an invalid gene index, so in wild
            # type experiment only genes that regulates no other genes are considered as output variables
            if regulates.issubset({perturbed_gene}):
                output_variables.add(act_gene)

        return output_variables.intersection(input_variables)

    def remove_perturbed_gene(self, perturbed_gene: int) -> bool:
        if perturbed_gene == -1:
            return False
        return len(self.get_regulated_by(perturbed_gene)) > 0
