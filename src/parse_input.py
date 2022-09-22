from typing import Dict, Optional, List, Set, Tuple
from classes import Regulation


State = Tuple[bool]


def read_input_matrix(matrix_path: str, tsv: bool) -> Dict[int, List[State]]:
    """Reads input matrix of steady-states, returns dictionary where key corresponds to perturbed gene or -1 for
    wild type experiment; and value corresponds to set of observed steady states of the target network.
    :param matrix_path - path to file containing the input matrix of steady states
    :param tsv - True if matrix's delimiter is tab, False if semicolon is used
    :return dictionary where key is index of perturbed gene and value is set of corresponding steady states
    TODO allow to have both - KO and OE experiments for the same gene, distinguish them in the input dictionary"""

    with open(matrix_path, "r") as matrix:
        steady_states: Dict[int, List[State]] = {}
        for line in matrix:
            perturbed_gene, state = parse_line(line, tsv)
            steady_states.setdefault(perturbed_gene, list()).append(state)
    return steady_states


def parse_line(line: str, tsv: bool) -> Tuple[int, State]:
    """Parses one input line and returns tuple containing index of perturbed gene as first element and corresponding
    state as second element.
    :param line - one input string line
    :param tsv - True if matrix's delimiter is tab, False if semicolon is used"""

    tmp = [int(x) for x in line.split(None if tsv else ";")]
    return tmp[0], tuple(bool(x) for x in tmp[1:])


def read_input_constrains(constraints_path: str, tsv: bool) -> Set[Regulation]:
    """Reads input file containing fixed regulations, returns list of initialized regulation objects
    :param constraints_path - path to file containing the fixed regulations' notations
    :param tsv - True if matrix's delimiter is tab, False if semicolon is used
    :return list of fixed regulations"""

    with open(constraints_path, "r") as constraints:
        regulations: Set[Regulation] = set()
        for line in constraints:
            tmp = tuple(line.split(None if tsv else ";"))
            source, target, sign = int(tmp[0]), int(tmp[1]), True if tmp[2] == '+' else False
            regulations.add(Regulation(source, target, sign))
    return regulations
