import src.classes.Generation as gen
import os


def output_to_directory(target_dir: str, gen: gen.Generation, num_of_vars: int, net_index: int):
    """Creates output directory from given path, creates output files and stores best networks each to separate file

    :param target_dir   output dir if not exists then will be created
    :param gen          last generation of networks
    :param num_of_vars  number of variables
    :param net_index    index of input network"""

    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)
    best_score = gen.best
    file_index = 0
    for i in range(gen.num_of_nets):
        if gen.scores[i] == best_score:
            with open("{}/{}_nodes_{}_instance_{}".format(target_dir, num_of_vars, net_index, file_index),
                      "w") as out_net_f:
                print(gen.networks[i].to_aeon_string(-1, set()), file=out_net_f)
            file_index += 1

    return best_score
