import numpy as np
np.random.seed(19680801)
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

sim = [('A', 0.999),
       ('B', 0.999),
       ('C', 0.999),
       ('H', 0.999),
       ('I', 0.999),
       ('J', 0.999),
       ('IJ', 0.999),
       ('D', 0.998),
       ('E', 0.998),
       ('F', 0.998),
       ('CD', 0.99),
       ('HIJ',0.99),
       ('AB',0.99),
       ('G',0.98),
       ('GHIJ',0.95),
       ('EF',0.94),
       ('EFGHIJ',0.94),
       ('ABCD',0.9)]

pangtree = [('A', 1),
            ('B', 1),
            ('C', 1),
            ('H', 1),
            ('I', 1),
            ('J', 1),
            # ('IJ', 0), #takiego wierzchołka nie ma
            ('D', 1),
            ('E', 1),
            ('F', 1),
            ('CD', 0.9954),
            ('HIJ',0.9953),
            ('AB',0.9968),
            ('G',1),
            ('GHIJ',0.9534),
            ('EF',0.9920),
            ('EFGHIJ',0.8420),
            ('ABCD',0.9507)]
x = [d[0] for d in sim]
y = [d[1] for d in sim]

x_pangtree = [d[0] for d in pangtree]
y_pangtree = [d[1] for d in pangtree]
ax.scatter(x, y, c="red", label="Evolution tree - distance to root")
ax.scatter(x_pangtree, y_pangtree, c="blue", label="Pangtree result - compatibility")

ax.legend(loc='lower left')
# ax.grid(True)
# plt.yticks(np.arange(0, max([max(y), max(y_pangtree)])+0.01, 0.01))
plt.xlabel("Node ID")
plt.xticks(rotation=90)
plt.title("Evolution tree vs. Pangtree result - comparison")
plt.savefig("tree_distances_tree3_nocycles500_p_1_stop_0999.png")
plt.show()
