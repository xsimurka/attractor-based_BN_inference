import os
from random import randint, random
from operator import itemgetter
from typing import Optional, List


def explicit_sinks_enumerator(model_path: str, output_path: str, num_of_vars: int,
                              enumerator_path: Optional[str] = ".\\sink-state-enumerator.exe",
                              reduction: Optional[float] = None) -> None:
    """The main function of this script. The script enumerates explicit sinks of given model using given enumerator
    to the output file specified.
    :param model_path       total file path to model that is desired to be processed
    :param output_path      total file path to desired output file with enumerated sinks
    :param num_of_vars      number of variables of given BN in <model_path>
    :param enumerator_path  total file path to sink state enumerator binary
    :param reduction        real number from [0;1] that determines the reduction of all explicit sinks
                            that are printed to <output_path> file. If None given, all sinks are printed"""

    with open(model_path, "r") as model_f:
        model_str = model_f.read()
    sinks = get_model_explicit_sinks(model_str, enumerator_path)
    write_sinks_to_file(sinks, output_path, -1, reduction)
    for i in range(num_of_vars):
        pert_model = get_perturbed_model(model_str, i, False)
        sinks = get_model_explicit_sinks(pert_model, enumerator_path)
        write_sinks_to_file(sinks, output_path, i, reduction)

        pert_model = get_perturbed_model(model_str, i, True)
        sinks = get_model_explicit_sinks(pert_model, enumerator_path)
        write_sinks_to_file(sinks, output_path, i, reduction)


def get_model_explicit_sinks(model: str, enumerator_path: Optional[str] = ".\\sink-state-enumerator.exe") -> List[List[int]]:
    """Functions returns all explicit sinks of the given model in .aeon format
    :param model            model in .aeon format
    :param enumerator_path  total file path to sink state enumerator binary
    :return                 list of model's explicit sinks"""

    r1 = randint(1000, 9999)
    tmp1 = "tmp_{}".format(r1)
    r2 = randint(100, 999)
    tmp2 = "tmp_{}".format(r2)
    with open(tmp2, "w") as model_f:
        print(model, file=model_f)
    command = "type {} | {} > {}".format(tmp2, enumerator_path, tmp1)
    os.system(command)
    os.remove(tmp2)

    result: List[List[int]] = []
    with open(tmp1, "r") as enum_out_f:
        for line in enum_out_f:
            if "Explicit sinks" in line:
                break

        for line in enum_out_f:
            # parses state from string representation and sort according to variable indices
            state: List[int] = list(map(lambda x: int(x[1] == 'true'),
                                        sorted(map(lambda x: (int(x[0][2:]), x[1]),
                                                   map(lambda x: x.split(": "), line.split("; ")[:-1])),
                                               key=itemgetter(0))))
            result.append(state)
    os.remove(tmp1)
    return result


def get_perturbed_model(model: str, pert_gene: int, pert_type: bool) -> str:
    """Function returns perturbed version of given model depending on given parameters
    :param model      model in .aeon format
    :param pert_gene  id of perturbed gene
    :param pert_type  True for over-expression, False for knockout
    return:           perturbed model in .aeon format"""

    result = str()
    lines = model.split('\n')
    for line in lines:
        if line:
            # regulation targeting perturbed gene is set to non-observable
            if "-| v_{}".format(pert_gene) in line or "-> v_{}".format(pert_gene) in line:
                result += replace_for_non_observable(line) + '\n'

            # update function of perturbed gene is skipped here...
            elif "$v_{}:".format(pert_gene) in line:
                pass

            else:
                result += line + '\n'

    # ...and replaced here for the perturbation type
    result += "$v_{}:{}\n".format(pert_gene, str(pert_type).lower())
    return result


def replace_for_non_observable(reg: str) -> str:
    """Function replaces the observable regulation for the unobservable.
    :param reg  observable regulation
    :return     unobservable variant"""

    s = reg.split(" -| ")
    if len(s) == 2:
        return s[0] + " -|? " + s[1]

    s = reg.split(" -> ")
    return s[0] + " ->? " + s[1]


def write_sinks_to_file(sinks: List[List[int]], output_path: str, pert_gene: int, reduction: float) -> None:
    """Function writes given sinks to given file
    :param sinks        list of explicit sinks
    :param output_path  total path to the output file
    :param pert_gene    perturbed gene id
    :param reduction    real number from [0;1] that determines the reduction of all explicit sinks
                        that are printed to <output_path> file. If None given, all sinks are printed"""

    with open(output_path, "a") as output_f:
        for state in sinks:
            if reduction is not None and random() > reduction:
                continue

            print(pert_gene, end='\t', file=output_f)
            for var in state:
                print(var, end='\t', file=output_f)
            print(file=output_f)
