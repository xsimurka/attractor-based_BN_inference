from typing import Set
import src.utils as utils
from src.classes.Regulation import Regulation
from parse_input import read_input_matrix, read_input_constrains
from output import output_to_directory
from src.classes.BNInfo import BNInfo


def main(sink_matrix_path: str, input_constraints_path: str, tsv: bool, threshold: float, best_ratio: float,
         num_of_nets: int, num_of_variables: int, num_of_genes: int, num_of_mutations: int, max_iter: int,
         max_fit: float, output_path: str):

    sinks: BNInfo = read_input_matrix(sink_matrix_path, tsv)  # DONE
    input_constraints: Set[Regulation] = read_input_constrains(input_constraints_path, tsv)  # DONE
    # derived_constraints: Set[Regulation] = derive_constraints(sinks, threshold)  # TODO
    act_generation = utils.create_initial_generation(num_of_nets, num_of_variables, input_constraints, set(), sinks)

    act_iter = 1
    while True:
        print("Actual iteration: ", act_iter)
        act_generation.compute_fitness()
        print(act_generation.best)
        if act_iter > max_iter or act_generation.best > max_fit:
            output_to_directory(output_path, act_generation)
            return

        act_generation = act_generation.create_new_generation(num_of_genes, num_of_mutations, best_ratio)
        act_iter += 1


if __name__ == '__main__':
    main("boolean_matrix_example", "input_constraints_example", True, 0.0, 0.2, 10, 4, 1, 1, 10, 0.95, "out_dir")
