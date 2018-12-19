from timeit import timeit

l = "".join([f'L{i}' for i in range(0, 10000, 4)])
s = "".join([f'S{i}' for i in range(0, 5000, 5)])
a = "A100"

line = l + s + a


def foo1(line):
    line = line[2:]
    line_iter = iter(line)

    node_parameters = {'L': [], 'S': [], 'A': [], '': []}
    current_node_parameter = next(line_iter)
    number_start = 1
    number_end = 1
    for i, c in enumerate(line_iter):
        if c == 'L' or c == 'S' or c == 'A':
            node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
            current_node_parameter = c
            number_start = i + 2
            number_end = i + 2
        else:
            number_end += 1
    node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
    aligned_to = node_parameters['A'][0] if node_parameters['A'] else None
    return node_parameters['L'], node_parameters['S'], aligned_to


# def foo2(line):

def foo_iter(line):
    all_positions = [i for i, s in enumerate(line) if s == 'L']

def foo_re(line):
    startingPos = 0

    position = line.find("L", startingPos)
    while position > -1:
        print(position)
        startingPos = position + 1
        position = text.find("|", startingPos)


# print(line)
# print(foo1(line))
print(f"foo1: {timeit('foo1(line)', globals=globals())}")
# print(f"foo2: {timeit('foo2(line)', globals=globals())}")
