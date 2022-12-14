from typing import List, Set, Tuple, Optional
import src.classes.Regulation as reg
import src.classes.TargetBN as bn

State = Tuple[bool]


def read_input_matrix(matrix_path: str, num_of_variables: int, delimiter: Optional[str], inputs, outputs) -> bn.TargetBN:
    """Reads input matrix of steady-states, returns dictionary where key corresponds to perturbed gene or -1 for
    wild type experiment; and value corresponds to set of observed steady states of the target network.

    :param matrix_path - path to file containing the input matrix of steady states
    :param num_of_variables
    :param delimiter - delimiter of matrix elements, if None given, whitespace characters are used
    :param inputs
    :param outputs
    :return dictionary where key is index of perturbed gene and value is set of corresponding steady states"""

    print("Reading input steady states...", end='')
    target_bn_info = bn.TargetBN(num_of_variables, inputs, outputs)
    with open(matrix_path, "r") as matrix:
        for line in matrix:
            perturbed_gene, state = parse_line(line, delimiter)
            get_experiment_state_list(target_bn_info, perturbed_gene, state).append(state)
    print(" done.")
    return target_bn_info


def parse_line(line: str, delimiter: str) -> Tuple[int, State]:
    """Parses one input line and returns tuple containing index of perturbed gene as first element and corresponding
    state as second element.

    :param line       one input string line
    :param delimiter  delimiter of matrix elements
    :return           tuple index of perturbed gene (or -1 for wild type), parsed Boolean vector"""

    tmp = [int(x) for x in line.split(delimiter)]
    return tmp[0], tuple(bool(x) for x in tmp[1:])


def read_input_constrains(constraints_path: str, delimiter: Optional[str]) -> Set[reg.Regulation]:
    """Reads input file containing fixed regulations, returns list of initialized regulation objects.

    :param constraints_path  path to file containing the fixed regulations' notations
    :param delimiter         delimiter of matrix elements
    :return                  list of fixed regulations"""

    print("Reading input constraints...", end="")
    with open(constraints_path, "r") as constraints:
        regulations: Set[reg.Regulation] = set()
        for line in constraints:
            tmp = tuple(line.split(delimiter))
            source, target, sign = int(tmp[0]), int(tmp[1]), True if tmp[2] == '+' else False
            regulations.add(reg.Regulation(source, target, sign))
    print(" done.")
    return regulations


def get_experiment_state_list(bn_info: bn.TargetBN, perturbed_gene: int, state: State) -> List[State]:
    """Returns list of steady-states of desired <perturbed_gene> according to given state, or an empty list if given
    gene has no steady-states of such experiment yet.
    If <state> has on <perturbed_gene> position 0 then list of knock-outs of given gene is returned (if exists)
    If <state> has on <perturbed_gene> position 1 then list of over-expressions of given gene is returned (if exists)
    If <perturbed_gene> is -1 then list of wild-types is returned

    :param bn_info         info about actual target network
    :param perturbed_gene  index of perturbed gene
    :param state           state that is currently being parsed from input
    :return                appropriate list of steady-states where the newly parsed steady-state should be added
                           if such list does not exist in corresponding dictionary, then new record is added"""

    if perturbed_gene == -1:
        return bn_info.wt_sinks
    experiment_dct = bn_info.oe_sinks if state[perturbed_gene] else bn_info.ko_sinks
    return experiment_dct.setdefault(perturbed_gene, list())
