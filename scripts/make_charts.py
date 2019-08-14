import sys
from pathlib import Path
import csv
from typing import List, Dict
import matplotlib.pyplot as plt
import numpy as np

data = Path(sys.argv[1]).resolve()
dataset_name = sys.argv[2]
plot_path = Path(sys.argv[3]).resolve()
data_output_path = Path(sys.argv[4]).resolve()


def make_plot():
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

    p_memory_time = []
    for p in p_to_memory_values.keys():
        p_memory_time.append((p, np.mean(p_to_memory_values[p]), np.mean(p_to_time_values[p])))
    # X Axis - Memory
    # Y Axis - Time
    ax.scatter([d[1] for d in p_memory_time], [d[2] for d in p_memory_time])
    for d in p_memory_time:
        ax.annotate(d[0], (d[1], d[2]))

    ax.set(xlabel='Memory [kB]', ylabel='Time [s]', title=f"Resources usage")
    fig.savefig(plot_path)

    with open(data_output_path, mode='w') as csv_file:
        resources_usage_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        resources_usage_writer.writerow(["P", "Memory", "Time"])
        for d in p_memory_time:
            resources_usage_writer.writerow([d[0], d[1], d[2]])



make_plot()