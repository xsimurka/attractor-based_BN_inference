from copy import deepcopy
from typing import Set, Optional, List
from random import choices
from src.classes.UpdateFunction import UpdateFunction
from src.classes.Regulation import Regulation
from src.detect import detect_steady_states, get_model_computational_structures, is_attractor_state
import src.utils as utils
import src.classes.BNInfo as bn
import src.classes.BipartiteGraph as bg


class BooleanNetwork:
    """Class represents Boolean network (BN)

    Attributes:
        target_bn_info    object that holds the information about target BN
        functions         list of <num_of_variables> update functions in the form of NCFs"""

    def __init__(self, target_bn_info: bn.BNInfo):
        self.target_bn_info = target_bn_info
        self.functions = [UpdateFunction(i, target_bn_info) for i in range(target_bn_info.num_of_vars)]

    def __deepcopy__(self):
        result: BooleanNetwork = BooleanNetwork(self.target_bn_info)
        result.functions = deepcopy(self.functions)
        return result

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
        :param fixed       True if regulation can not be removed afterwards, otherwise False
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

        non_input_genes = list(set(range(self.target_bn_info.num_of_vars)) - self.target_bn_info.input_genes)
        genes_to_mutate = choices(non_input_genes, k=num_of_genes)
        for gene in genes_to_mutate:
            self.functions[gene].mutate(num_of_mutations, self.total_num_of_regulations)

    def to_aeon_string(self, perturbed_gene: int, isolated_variables: Set[int],
                       perturbation_state: Optional[bool] = None) -> str:
        """Returns string representation of given BN in .aeon format. Dimension of the model can be reduced due
        to the existence of isolated variables that will not be included in attractor analysis.

        :param perturbed_gene      denotes perturbed gene index, None if there is non perturbation
        :param isolated_variables  set of variables that do not regulate other variables and are not regulated too
        :param perturbation_state  False if perturbation is a knockout (KO) True if over-expression (OE),
                                   None if there is no perturbation
        :return                    AEON string representation of reduced model"""

        model_string = str()

        # adds all meaningful regulations to the model
        for i in range(self.target_bn_info.num_of_vars):  # iterates over update functions of given BN
            # perturbed gene acts as fixed input node therefore, its original regulators are skipped
            if i == perturbed_gene:
                continue

            for j in range(self.functions[i].arity):  # iterates over variables of i-th update function
                reg_type = ">" if self.functions[i].canalyzing[j] == self.functions[i].canalyzed[j] else "|"
                model_string += "v_{} -{} v_{}\n".format(self.functions[i].regulators[j], reg_type, i)

        # adds all meaningful update functions to the model
        for i in range(self.target_bn_info.num_of_vars):
            if i in isolated_variables:  # isolated variables are not added at all
                continue

            if i == perturbed_gene:  # perturbed gene has fixed update function - true (KO), false (OE)

                if self.functions[i].regulators != [i]:
                    # edge case: fixing variable that has only reflexive regulation results in its isolation
                    model_string += "$v_{}:{}\n".format(i, str(perturbation_state).lower())

            else:  # create update function iff gene is not perturbed and is not isolated
                model_string += "{}\n".format(self.functions[i].to_aeon_string())

        return model_string

    def compute_fitness(self, perturbed_gene_index: int,
                        perturbed_gene_state: Optional[bool] = None) -> float:
        """Evaluates fitness of given BN depending on its steady-state attractors comparing to steady-states from
        given attractor data.

        :param perturbed_gene_index   index of perturbed gene of currently evaluated BN (-1 if wild-type)
        :param perturbed_gene_state   True if current experiment is over-expression of <perturbed_gene_index>
                                      False if current experiment is knock-out of <perturbed_gene_index>
                                      None if current experiment is wild-type
        :return                       real number from [0;1] that determines BN's fitness to input data"""

        isolated_variables = self.get_isolated_variables(perturbed_gene_index)
        aeon_model_string = self.to_aeon_string(perturbed_gene_index, isolated_variables, perturbed_gene_state)
        model, sag = get_model_computational_structures(aeon_model_string)
        target_sinks = self.target_bn_info.get_target_sinks(perturbed_gene_index, perturbed_gene_state)

        if isolated_variables:  # reduce dimension of target sinks in order to match observed data dimension
            target_sinks = utils.reduce_attractors_dimension(target_sinks, isolated_variables)

        observed_sinks = detect_steady_states(sag)
        bpg = bg.BipartiteGraph(target_sinks, observed_sinks)
        cost, pairs = bpg.minimal_weighted_assignment()
        dimension = self.target_bn_info.num_of_vars - len(isolated_variables)  # number of variables in each state after reduction
        matching_variables = 0  # variable = one element of a state vector

        # if some states were matched, then total number of matching variables is equal to total number of variable
        # pairs minus cost (total number of variables that do not match), else it stays equal to 0 (initial value)
        if cost is not None:
            matching_variables += (min(len(target_sinks), len(observed_sinks)) * dimension) - cost
        total_variables = len(target_sinks) * dimension

        if len(target_sinks) > len(observed_sinks):
            # some steady-states from data were not reached by actual model, which is problematic because such
            # model can not capture desired behaviour specified by the input data
            # total_variables is therefore, equal to number of variable pairs plus number of unpaired variables

            unmatched_states = utils.get_unmatched_states(bpg, pairs, target_sinks, position=0)
            # check if missing steady-states are on some other type of attractor because there is possibility of
            # breaking such attractor to steady-state by tiny mutation of actual BN
            for state in unmatched_states:
                if is_attractor_state(model, sag, state):
                    matching_variables += dimension * 1 / 2  # penalty 1/2 for not being in single state attractor

        # else:  len(target_sinks) <= len(observed_sinks)
            # there is possibility that some steady-states were not caught while measuring which is more and more
            # probable by increasing number of genes where the number of steady-states grows exponentially
            # therefore, no penalty is given and total number of variables is only equal to number of variable pairs


        return matching_variables / total_variables

    def get_regulated_by(self, gene: int) -> Set[int]:
        """Returns set of genes that are regulated by given gene.

        :param gene  target regulator
        :return      set of genes regulated by given gene"""

        result = set()
        for i in range(self.target_bn_info.num_of_vars):
            for j in self.functions[i].regulators:
                if gene == j:
                    result.add(i)

        return result

    def get_isolated_variables(self, perturbed_gene: int) -> Set[int]:
        """Returns set of isolated variables of given BN i.e., variables that have no regulators (inputs)
         and do not regulate any other variable (outputs) at the same time.

        :return  set of isolated variables"""

        # perturbed gene is definitely input because it looses its update function so also its regulators
        input_variables = set() if perturbed_gene == -1 else {perturbed_gene}
        output_variables = set()

        # input variables
        for act_gene in range(self.target_bn_info.num_of_vars):
            if not self.functions[act_gene].regulators:  # only genes that have empty update function are inputs
                input_variables.add(act_gene)

        # output variables
        for act_gene in range(self.target_bn_info.num_of_vars):
            regulates = self.get_regulated_by(act_gene)

            # gene is an output variable iff it regulates either no other genes or only perturbed gene (this regulation
            # will be formally removed in mutant type experiment by fixing expression level of perturbed gene)
            # in wild type experiment, <perturbed_gene> is equal to -1, which is an invalid gene index, so in wild
            # type experiment only genes that regulates no other genes are considered as output variables
            if regulates.issubset({perturbed_gene}):
                output_variables.add(act_gene)

        return output_variables.intersection(input_variables)
