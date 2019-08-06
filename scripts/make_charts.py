import sys
from pathlib import Path
import csv
from typing import List, Dict
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

data = Path(sys.argv[1]).resolve()
dataset_name = sys.argv[2]
plot_path = Path(sys.argv[3]).resolve()
data_output_path = Path(sys.argv[4]).resolve()

def version1(memory_plot_path, time_plot_path):
    p_to_memory_values: Dict[float, List[int]] = {}
    p_to_time_values: Dict[float, List[int]] = {}

    with open(data) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader, None)
        line_count = 0
        for row in csv_reader:
            p = row[2]
            if row[2] in p_to_memory_values:
                p_to_memory_values[p].append(int(row[4]))
                p_to_time_values[p].append(float(row[3]))
            else:
                p_to_memory_values[p] = [int(row[4])]
                p_to_time_values[p] = [float(row[3])]

    fig, ax = plt.subplots()
    for p, memory_values in p_to_memory_values.items():
        ax.plot(list(range(len(memory_values))), memory_values, label=p)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set(xlabel='Iteration', ylabel='Memory [kB]', title=f"Memory usage - dataset {dataset_name}")
    ax.legend(loc=1)
    ax.grid()
    plt.xticks(list(range(len(p_to_memory_values["1"]))))
    fig.savefig(memory_plot_path)

    fig, ax = plt.subplots()
    for p, time_values in p_to_time_values.items():
        ax.plot(list(range(len(time_values))), time_values, label=p)

    ax.set(xlabel='Iteration', ylabel='Time [s]', title=f"Time usage - dataset {dataset_name}")
    plt.xticks(list(range(len(p_to_time_values["1"]))))
    ax.legend(loc=1)
    ax.grid()
    fig.savefig(time_plot_path)


def version2():
    p_to_memory_values: Dict[float, List[int]] = {}
    p_to_time_values: Dict[float, List[int]] = {}

    with open(data) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader, None)
        for row in csv_reader:
            p = row[2]
            if row[2] in p_to_memory_values:
                p_to_memory_values[p].append(int(row[4]))
                p_to_time_values[p].append(float(row[3]))
            else:
                p_to_memory_values[p] = [int(row[4])]
                p_to_time_values[p] = [float(row[3])]

    fig, ax = plt.subplots()

    p_m_t = []
    for p in p_to_memory_values.keys():
        p_m_t.append((p, np.mean(p_to_memory_values[p]), np.mean(p_to_time_values[p])))
    # X Axis - Memory
    # Y Axis - Time
    ax.scatter([d[1] for d in p_m_t], [d[2] for d in p_m_t])
    for d in p_m_t:
        ax.annotate(d[0], (d[1], d[2]))

    ax.set(xlabel='Memory [kB]', ylabel='Time [s]', title=f"Resources usage")
    fig.savefig(plot_path)

    with open(data_output_path, mode='w') as csv_file:
        resources_usage_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        resources_usage_writer.writerow(["P", "Memory", "Time"])
        for d in p_m_t:
            resources_usage_writer.writerow([d[0], d[1], d[2]])



version2()