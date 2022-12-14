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
    for i in range(gen.num_of_nets):
        if gen.scores[i] == best_score:
            with open("{}\\model{}_{}.aeon".format(target_dir, num_of_vars, net_index),
                      "w") as out_net_f:
                aeon_str = gen.networks[i].to_aeon_string(-1, set())
                print(aeon_str, file=out_net_f)
