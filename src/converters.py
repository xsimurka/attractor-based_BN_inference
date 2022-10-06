from re import sub
from typing import Optional


# import src.classes.BooleanNetwork as bn


def zavorky(s):
    result = 0
    for char in s:
        if char == '(':
            result += 1
        if char == ')':
            result -= 1
    return result == 0


def file_converter(input_file: str, output_format: str):
    if output_format.lower() == "boolesim":
        output_file = sub(r'^(.+)\.([^.]+)$', r'\1_bsim.\2', input_file)
        if output_file == input_file:
            output_file += "_output.txt"
        to_boolesim(input_file, output_file)
    elif output_format.lower() == "aeon":
        output_file = sub(r'^(.+)\.([^.]+)$', r'\1_aeon.\2', input_file)
        if output_file == input_file:
            output_file += "_output.aeon"
        to_aeon(input_file, output_file)
    else:
        print("Invalid output file format. Choose from {\"boolesim\", \"aeon\"}.")


def to_boolesim(input_path, output_path):
    with open(input_path, "r") as input_file, open(output_path, "w") as output_file:
        for line in input_file:
            line = sub('AND|and', '&&', line)
            line = sub('OR|or', '||', line)
            line = sub('NOT|not', '!', line)

            print(line, file=output_file, end='')


def to_aeon(input_path, output_path):
    with open(input_path, "r") as input_file, open(output_path, "w") as output_file:
        for line in input_file:
            line = sub('AND|and', '&', line)
            line = sub('OR|or', '|', line)
            line = sub('NOT|not', '!', line)

            line = sub(r"\s*=\s*", ":", line)
            i = len(line) - 1
            base = get_num(line, 0)

            lst = []
            while True:
                if line[i] == ':':
                    break

                if line[i].isdigit() and not line[i - 1].isdigit():
                    actual = get_num(line, i)
                    bond = ">" if line[i - 1] == ' ' else "|"
                    lst.append("v_{0} -{1} v_{2}".format(actual, bond, base))
                    line = line[:i] + "v_" + line[i:]
                i -= 1

            if line[-1] == ':':
                print("#position:v_{0}:0,0".format(line[:-1]), file=output_file)
            else:
                line = delete_redundant_brackets(line)
                line = refine_spaces(line)
                print("$v_" + line, file=output_file, end='' if line[-1] == '\n' else '\n')

            for elem in lst:
                print(elem, file=output_file)


def get_num(string: str, index: int) -> str:
    result = ''
    while index < len(string):
        if string[index].isdigit():
            result += string[index]
        else:
            break
        index += 1
    return result


def delete_redundant_brackets(line: str) -> str:
    old: Optional[str] = None
    while old != line:
        old = line
        line = sub(r"\(([^(|&]*)\)", r"\1", line)
    return line


def refine_spaces(line: str) -> str:
    line = sub(r'\s*([|&()])\s*', r" \1 ", line)
    line = sub(r'\s+', ' ', line)
    return line


def from_ncfs_to_aeon(bn) -> str:
    result = ""
    for i in range(len(bn.functions)):
        if not bn.functions[i].indices:  # empty update function
            result += "#position:v_{0}:0,0\n".format(i)
            continue
        for j in range(len(bn.functions[i].indices)):
            bond = ">" if bn.functions[i].canalyzing[j] == bn.functions[i].canalyzed[j] else "|"
            result += "v_{0} -{1} v_{2}\n".format(bn.functions[i].indices[j], bond, i)

    for i in range(len(bn.functions)):
        result += "{}\n".format(bn.functions[i])

    return result


file_converter("test_network", "aeon")
