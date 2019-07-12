import numpy as np
np.random.seed(19680801)
import matplotlib.pyplot as plt


fig, ax = plt.subplots()


data = [
    # (1 - 0.9968, 0.001),  # A - AB
    (1 - 1,      0.001),  # B - AB

    # (1 - 0.9954, 0.001),  # C - CD
    (1 - 1,      0.002),  # D - CD

    # (1 - 0.9987, 0.021),  # A - ABCD
    (1 - 0.9991, 0.021),  # B - ABCD
    # (1 - 0.9517, 0.011),  # C - ABCD
    # (1 - 0.9507, 0.012),  # D - ABCD

    # (1 - 0.992, 0.002),  # E - EF
    (1 - 1, 0.002),  # F - EF

    # (1 - 0.9953, 0.001),  # H - HIJ
    (1 - 1, 0.002),  # I - HIJ
    # (1 - 0.9966, 0.002),  # J - HIJ

    (1 - 0.9534, 0.02),  # G - GHIJ
    # (1 - 0.9966, 0.011),  # H - GHIJ
    # (1 - 0.999, 0.012),  # I - GHIJ
    # (1 - 0.9979, 0.012),  # J - GHIJ

    # (1 - 0.842, 0.062),  # E - EFGHIJ
    # (1 - 0.8424, 0.062),  # F - EFGHIJ
    (1 - 0.9534, 0.07),  # G - EFGHIJ
    # (1 - 0.9973, 0.061),  # H - EFGHIJ
    # (1 - 0.9983, 0.062), # I - EFGHIJ
    # (1 - 0.9979, 0.062),  # J - EFGHIJ

    # (1 - 0.7247, 0.121), # A - root
    # (1 - 0.7246, 0.121), # B - root
    # (1 - 0.7307, 0.111), # C - root
    # (1 - 0.7301, 0.112), # D - root
    # (1 - 0.842, 0.122), # E - root
    # (1 - 0.8424, 0.122), # F - root
    (1 - 0.9534, 0.13), # G - root
    # (1 - 0.9973, 0.121), # H - root
    # (1 - 0.9983, 0.122), # I - root
    # (1 - 0.9979, 0.122), # J - root

    (0, 0),  # G - G
]
# (((A:0.001,B:0.001)AB:0.02,(C:0.001,D:0.002)CD:0.01)
# ABCD:0.1,((E:0.002,F:0.002)EF:0.06,(G:0.02,(H:0.001,(I:0.001,J:0.001)IJ:0.001)HIJ:0.01)GHIJ:0.05)
# EFGHIJ:0.06);

x = [d[0] for d in data]
y = [d[1] for d in data]

ax.scatter(x, y, c="red")
plt.xlabel("1 - minComp")
plt.ylabel("Evolution tree distance")
plt.title("Evolution tree vs. Pangtree result - comparison")
plt.savefig("tree_distances_version3_tree3_nocycles500_p_1_SINGLENODES_smaller_x.png")
plt.show()
