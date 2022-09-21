from classes import Generation
import os


def output_to_directory(target_dir: str, gen: Generation):
    """"""

    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)
    best_score = gen.best
    file_index = 0
    for i in range(gen.num_of_nets):
        if gen.scores[i] == best_score:
            with open("output_network_{}".format(file_index), "w") as out_net_f:
                print(gen.networks[i].to_aeon_string(-1, set()), file=out_net_f)
            file_index += 1

    return best_score
