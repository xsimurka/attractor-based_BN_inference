from typing import Set
import src.utils as utils
from src.classes.Regulation import Regulation
from parse_input import read_input_matrix, read_input_constrains
from output import output_to_directory
from src.classes.BNInfo import BNInfo

# True if files above use tab as separator, otherwise False
TSV = True
# threshold for correlation constraints
THRESHOLD = 0.85
# how many % of best fitting BNs are automatically picked to the next generation
ELITE_RATIO = 0.2
# number of networks in each generation
NUM_OF_NETWORKS = 10
# number of genes that are mutated in each generation
NUM_OF_GENES = 2
# number of mutations performed on each mutated gene
NUM_OF_MUTATIONS = 2
# maximum number of iterations of genetic algorithm
MAX_ITERATION = 30
# when some network reach this fitness, algorithm ends immediately
MAX_FITNESS = 0.95


def main(sink_matrix_path: str, input_constraints_path: str, num_of_variables: int,
         output_path: str, inputs: Set[int], outputs: Set[int], net_index: int, output: bool):
    """Entry point of the whole inference algorithm.

    :param sink_matrix_path        path to the file with steady-state matrix
    :param input_constraints_path  path to the file with constraints
    :param num_of_variables        number of variables of each network
    :param output_path             path for directory containing output files
    :param inputs                  set of input nodes i.e., nodes that are not regulated by other genes
    :param outputs                 set of output nodes i.e., nodes that do not regulate other genes
    :param net_index               index of input net
    :param output                  True if output results to files, otherwise False"""

    print("Start")
    target_bn_info: BNInfo = read_input_matrix(sink_matrix_path, num_of_variables, TSV, inputs, outputs)
    input_constraints: Set[Regulation] = read_input_constrains(input_constraints_path, TSV)
    derived_constraints: Set[Regulation] = target_bn_info.derive_constraints(THRESHOLD)
    act_generation = utils.create_initial_generation(NUM_OF_NETWORKS, input_constraints,
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
        if act_iter >= MAX_ITERATION or act_generation.best >= MAX_FITNESS:
            print("Ending genetic algorithm...")
            print("Writing best fitting networks to files...", end="")
            if output:
                output_to_directory(output_path, act_generation, num_of_variables, net_index)
                print(" done.")
                return
            else:
                print(" done.")
                return act_generation


        act_generation = act_generation.create_new_generation(NUM_OF_GENES, NUM_OF_MUTATIONS, ELITE_RATIO)
        act_iter += 1


if __name__ == '__main__':
    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_1_constraints",
         10, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir", {1, 8, 9}, {0, 3, 5}, 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_constraints",
         10, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir", {2, 5, 9}, {1, 7, 8}, 2)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_1_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_1_constraints",
         20, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir", {10, 13, 15},
         {7, 9, 11, 12, 17, 18}, 1)

    main(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_2_sinks",
         r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\20_nodes_2_constraints",
         20, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir", set(), set(), 2)
