from typing import Dict, List, Set, Tuple
import src.classes.Regulation as reg


State = Tuple[bool]


class BNInfo:

    def __init__(self):
        self.wt_sinks: List[State] = []
        self.ko_sinks: Dict[int, List[State]] = {}
        self.oe_sinks: Dict[int, List[State]] = {}


def read_input_matrix(matrix_path: str, tsv: bool) -> BNInfo:
    """Reads input matrix of steady-states, returns dictionary where key corresponds to perturbed gene or -1 for
    wild type experiment; and value corresponds to set of observed steady states of the target network.
    :param matrix_path - path to file containing the input matrix of steady states
    :param tsv - True if matrix's delimiter is tab, False if semicolon is used
    :return dictionary where key is index of perturbed gene and value is set of corresponding steady states
    TODO allow to have both - KO and OE experiments for the same gene, distinguish them in the input dictionary"""

    target_bn_info = BNInfo()
    with open(matrix_path, "r") as matrix:
        for line in matrix:
            perturbed_gene, state = parse_line(line, tsv)
            get_experiment_state_list(target_bn_info, perturbed_gene, state).append(state)
    return target_bn_info


def parse_line(line: str, tsv: bool) -> Tuple[int, State]:
    """Parses one input line and returns tuple containing index of perturbed gene as first element and corresponding
    state as second element.

    :param line  one input string line
    :param tsv  True if matrix's delimiter is tab, False if semicolon is used
    :return     tuple index of perturbed gene (or -1 for wild type), parsed Boolean vector"""

    tmp = [int(x) for x in line.split(None if tsv else ";")]
    return tmp[0], tuple(bool(x) for x in tmp[1:])


def read_input_constrains(constraints_path: str, tsv: bool) -> Set[reg.Regulation]:
    """Reads input file containing fixed regulations, returns list of initialized regulation objects.

    :param constraints_path  path to file containing the fixed regulations' notations
    :param tsv               True if matrix's delimiter is tab, False if semicolon is used
    :return                  list of fixed regulations"""

    with open(constraints_path, "r") as constraints:
        regulations: Set[reg.Regulation] = set()
        for line in constraints:
            tmp = tuple(line.split(None if tsv else ";"))
            source, target, sign = int(tmp[0]), int(tmp[1]), True if tmp[2] == '+' else False
            regulations.add(reg.Regulation(source, target, sign))
    return regulations


def get_experiment_state_list(bn_info: BNInfo, perturbed_gene: int, state: State) -> List[State]:
    """"""
    if perturbed_gene == -1:
        return bn_info.wt_sinks
    experiment_dct = bn_info.oe_sinks if state[perturbed_gene] else bn_info.ko_sinks
    return experiment_dct.setdefault(perturbed_gene, list())
