from typing import List
import src.utils as utils
from src.classes.Regulation import Regulation
from parse_input import read_input_matrix, read_input_constrains
from output import output_to_directory
from src.classes.BNInfo import BNInfo


def main(sink_matrix_path: str, input_constraints_path: str, tsv: bool, threshold: float, best_ratio: float,
         num_of_nets: int, num_of_variables: int, num_of_genes: int, num_of_mutations: int, max_iter: int,
         max_fit: float, output_path: str):
    """Entry point of the whole inference algorithm.

    :param sink_matrix_path        path to the file with steady-state matrix
    :param input_constraints_path  path to the file with constraints
    :param tsv                     True if files above use tab as separator, otherwise False
    :param threshold               threshold for correlation constraints
    :param best_ratio              how many % of best fitting BNs are automatically picked to the next generation
    :param num_of_nets             number of networks in each generation
    :param num_of_variables        number of variables of each network
    :param num_of_genes            number of genes that are mutated in each generation
    :param num_of_mutations        number of mutations performed on each mutated gene
    :param max_iter                maximum number of iterations of genetic algorithm
    :param max_fit                 when some network reach this fitness, algorithm ends immediately
    :param output_path             path for directory containing output files"""

    print("Start")
    sinks: BNInfo = read_input_matrix(sink_matrix_path, num_of_variables, tsv)
    input_constraints: List[Regulation] = read_input_constrains(input_constraints_path, tsv)
    derived_constraints: List[Regulation] = sinks.derive_constraints(threshold)
    act_generation = utils.create_initial_generation(num_of_nets, num_of_variables, input_constraints,
                                                     derived_constraints, sinks)

    print("Starting genetic algorithm...")
    act_iter = 1
    while True:
        print("Actual iteration: ", act_iter)
        print("Evaluating fitness...", end="")
        act_generation.compute_fitness()
        print(" done.")
        print("Best fitness of generation: {}".format(act_generation.best))
        print()
        if act_iter > max_iter or act_generation.best >= max_fit:
            print("Ending genetic algorithm...")
            print("Writing best fitting networks to files...", end="")
            output_to_directory(output_path, act_generation)
            print(" done.")
            return

        act_generation = act_generation.create_new_generation(num_of_genes, num_of_mutations, best_ratio)
        act_iter += 1


if __name__ == '__main__':
    main("boolean_matrix_example", "input_constraints_example", True, 0.0, 0.2, 10, 4, 1, 1, 10, 0.95, "out_dir")
