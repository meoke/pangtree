import numpy as np
np.random.seed(19680801)
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

sim = [('A', 0.001),
       ('B', 0.001),
       ('C', 0.001),
       ('H', 0.001),
       ('I', 0.001),
       ('J', 0.001),
       ('IJ', 0.001),
       ('D', 0.002),
       ('E', 0.002),
       ('F', 0.002),
       ('CD', 0.01),
       ('HIJ',0.01),
       ('AB',0.02),
       ('G',0.02),
       ('GHIJ',0.05),
       ('EF',0.06),
       ('EFGHIJ',0.06),
       ('ABCD',0.1)]

pangtree = [('A', 0.0032),
            ('B', 0.0032),
            ('C', 0.0046),
            ('H', 0.0047),
            ('I', 0.0047),
            ('J', 0.0047),
            ('IJ', 0), #takiego wierzchołka nie ma
            ('D', 0.0046),
            ('E', 0.008),
            ('F', 0.008),
            ('CD', 0.0447),
            ('HIJ',0.0419),
            ('AB',0.0461),
            ('G',0.0466),
            ('GHIJ',0.1114),
            ('EF',0.15),
            ('EFGHIJ',0.1174),
            ('ABCD',0.2261)]
x = [d[0] for d in sim]
y = [d[1] for d in sim]

x_pangtree = [d[0] for d in pangtree]
y_pangtree = [d[1] for d in pangtree]
ax.scatter(x, y, c="red", label="Simulated")
ax.scatter(x_pangtree, y_pangtree, c="blue", label="Pangtree")

ax.legend(loc='upper left')
# ax.grid(True)
plt.yticks(np.arange(0, max([max(y), max(y_pangtree)])+0.01, 0.01))
plt.xticks(rotation=90)
plt.title("Distance between node and parent")
plt.savefig("distance_node_parent_tree3_p_1_stop_0999.png")
plt.show()
