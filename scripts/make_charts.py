import sys
from pathlib import Path
import csv
from typing import List, Dict
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

data = Path(sys.argv[1]).resolve()
dataset_name = sys.argv[2]
time_plot_path = Path(sys.argv[3]).resolve()
memory_plot_path = Path(sys.argv[4]).resolve()

p_to_memory_values: Dict[float, List[int]] = {}
p_to_time_values: Dict[float, List[int]] = {}

with open(data) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    next(csv_reader, None)
    line_count = 0
    for row in csv_reader:
        if row[2] in p_to_memory_values:
            p_to_memory_values[row[2]].append(int(row[4]))
            p_to_time_values[row[2]].append(float(row[3]))
        else:
            p_to_memory_values[row[2]] = [int(row[4])]
            p_to_time_values[row[2]] = [float(row[3])]

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
