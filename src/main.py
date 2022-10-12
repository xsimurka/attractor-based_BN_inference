from typing import Set
import src.utils as utils
from src.classes.Regulation import Regulation
from parse_input import read_input_matrix, read_input_constrains
from output import output_to_directory
from src.classes.BNInfo import BNInfo


def main(sink_matrix_path: str, input_constraints_path: str, tsv: bool, threshold: float, best_ratio: float,
         num_of_nets: int, num_of_variables: int, num_of_genes: int, num_of_mutations: int, max_iter: int,
         max_fit: float, output_path: str, inputs: Set[int], outputs: Set[int], net_index: int):
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
    :param output_path             path for directory containing output files
    :param inputs                  set of input nodes i.e., nodes that are not regulated by other genes
    :param outputs                 set of output nodes i.e., nodes that do not regulate other genes
    :param net_index               index of input net """

    print("Start")
    target_bn_info: BNInfo = read_input_matrix(sink_matrix_path, num_of_variables, tsv, inputs, outputs)
    input_constraints: Set[Regulation] = read_input_constrains(input_constraints_path, tsv)
    derived_constraints: Set[Regulation] = target_bn_info.derive_constraints(threshold)
    act_generation = utils.create_initial_generation(num_of_nets, input_constraints,
                                                     derived_constraints - input_constraints, target_bn_info)

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
            output_to_directory(output_path, act_generation, num_of_variables, net_index)
            print(" done.")
            return

        act_generation = act_generation.create_new_generation(num_of_genes, num_of_mutations, best_ratio)
        act_iter += 1


if __name__ == '__main__':
    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_1_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 2)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_1_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_2_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 2)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\30_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\30_nodes_1_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\30_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\30_nodes_2_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 2)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\50_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\50_nodes_1_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\50_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\50_nodes_2_constraints",
         True, 0.85, 0.2, 10, 4, 1, 1, 10, 0.95, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir",
         set(), set(), 2)
