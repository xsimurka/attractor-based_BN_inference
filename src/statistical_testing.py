from main import *
import matplotlib.pyplot as plt
import pandas as pd
import classes.BooleanNetwork as bn

REPS = 100


def statistics_testing(sink_matrix_path: str, input_constraints_path: str, num_of_variables: int,
                       output_path: str, inputs: Set[int], outputs: Set[int], net_index: int):
    reg_counter = {}
    for i in range(num_of_variables):
        for j in range(num_of_variables):
            reg_counter[str(i) + "->" + str(j)] = 0
            reg_counter[str(i) + "-|" + str(j)] = 0

    for rep in range(REPS):
        print("Rep ", rep)
        g = main(sink_matrix_path, input_constraints_path, num_of_variables, output_path, inputs, outputs, net_index,
                 False)
        scores = pd.Series(g.scores)
        best = scores[scores == g.best].index.values.tolist()
        for index in best:
            count_regulations(g.networks[index], reg_counter)
    nbest = pd.Series(list(reg_counter.values())).nlargest(30).index.values.tolist()
    nkeys = list(reg_counter.keys())
    nvalues = list(reg_counter.values())
    courses = [nkeys[i] for i in nbest]
    values = [nvalues[i] for i in nbest]

    plt.figure(figsize=(50, 50))

    # creating the bar plot
    plt.bar(courses, values, width=0.4)

    plt.xlabel("Regulations")
    plt.ylabel("Number of occurrences")
    plt.title("10_nodes_1 100 reps statistics")
    plt.show()


def count_regulations(network: bn.BooleanNetwork, counter):
    for k in range(len(network.functions)):
        for j in network.functions[k].regulators:
            i, o = network.functions[k].get_c_and_c_values(j)
            sign = ">" if i == o else "|"
            counter["{}-{}{}".format(j, sign, k)] += 1


statistics_testing(r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_sinks",
                   r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\input_networks\10_nodes_2_constraints",
                   10, r"C:\Users\Andrej Šimurka\Desktop\atractor_analysis_inference\out_dir", {2, 5, 9}, {1, 7, 8}, 2)
