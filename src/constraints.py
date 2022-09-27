from math import sqrt
from typing import List, Tuple, Optional, Dict, Set
from src.classes.Regulation import Regulation


State = Tuple[bool]


def read_matrix(f: str):
    with open(f, "r") as matrix_file:
        result = []
        next(matrix_file)
        for line in matrix_file:
            result.append(list(map(lambda x: float(x), line.split()[1:])))

    return result


def read_structure(f: str):
    with open(f, "r") as struct_file:
        pos, neg = [], []
        for line in struct_file:
            source, target, bond = int(line.split()[0][1:]), \
                                   int(line.split()[1][1:]), \
                                   True if line.split()[2] == '+' else False
            if bond:
                pos.append((source, target))
            else:
                neg.append((source, target))
    return pos, neg


def count_diffs(all, bonds, matrix):
    bonds = set(bonds)
    max_not = -1
    max_not_s, max_not_t = -1, -1
    for source, target in all:
        s_wt = matrix[0][source - 1]
        t_wt = matrix[0][target - 1]
        s_ex = matrix[source][source - 1]
        t_ex = matrix[source][target - 1]
        value = abs(t_ex - t_wt) / abs(s_ex - s_wt)
        if (source, target) not in bonds and source != target and s_wt > 0.1:
            if value > max_not:
                max_not = value
                max_not_s = source
                max_not_t = target

        print('* ' if (source, target) in bonds else '  ', source, " -> ", target, ": ", value, sep='')
    print()
    print(max_not_s, " -> ", max_not_t, ": ", max_not, sep='')


def pearson_correlation(matrix):
    num_of_genes = len(matrix[0])
    result = [[1.0 for _ in range(num_of_genes)] for _ in range(num_of_genes)]
    average_expr = []
    for i in range(num_of_genes):
        sum = 0
        for exp in range(len(matrix)):
            if i != exp - 1:
                sum += matrix[exp][i]
        average_expr.append(sum/(len(matrix) - 1))

    for i in range(num_of_genes):
        for j in range(i + 1, num_of_genes):
            result[i][j] = result[j][i] = calc_corr(i, j, matrix, average_expr)

    return result


def calc_corr(k, l, matrix, average_expr):
    s1 = s2 = s3 = 0
    for exp in matrix:
        s1 += (exp[k] - matrix[0][k]) * (exp[l] - matrix[0][l])
        s2 += (exp[k] - matrix[0][k]) ** 2
        s3 += (exp[l] - matrix[0][l]) ** 2

    return s1 / (sqrt(s2) * sqrt(s3))


def constraints(corr_matrix: List[List[float]], threshold: float) -> Set[Regulation]:
    result: Set[Regulation] = set()
    for i in range(len(corr_matrix)):
        for j in range(len(corr_matrix[i])):
            if i == j:
                continue
            if abs(corr_matrix[i][j]) > threshold:
                result.add(Regulation(i, j, True if corr_matrix[i][j] > 0 else False))

    return result


def derive_constraints(sinks: Dict[int, List[State]], threshold: float) -> Set[Regulation]:
    # vezme slovnik spravi maticu a zavola 'constraints()'
    return set()


matrix = read_matrix("InSilicoSize10-Ecoli1-nonoise-null-mutants.tsv")
corr = pearson_correlation(matrix)
pos, neg = read_structure("InSilicoSize10-Ecoli1.tsv")

gen = [(i, j) for i in range(1,11) for j in range(1,11)]
count_diffs(gen, pos + neg, matrix)
