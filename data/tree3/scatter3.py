import numpy as np
np.random.seed(19680801)
import matplotlib.pyplot as plt


fig, ax = plt.subplots()


data = [
    (0, 0.121),  # A - A
    # (1-0.9968, 0.12),  # A - AB
    # (1-0.9987, 0.1),  # A - ABCD

    (0, 0.121),  # B - B
    (1 - 1, 0.12),  # B - AB
    (1 - 0.9991, 0.1),  # B - ABCD

    (0, 0.001),  # C - C
    # (1 - 0.9954, 0.11),  # C - CD
    # (1 - 0.9517, 0.1),  # C - ABCD

    (0, 0.002),  # D - D
    (1 - 1, 0.11),  # D - CD
    # (1 - 0.9507, 0.1),  # D - ABCD

    (0, 0.122),  # E - E
    # (1 - 0.992, 0.12),  # E - EF
    # (1 - 0.8420, 0.06),  # E - EFGHIJ

    (0, 0.122),  # F - F
    (1 - 1, 0.12),  # F - EF
    # (1 - 0.8424, 0.06),  # F - EFGHIJ

    (0, 0.13),  # G - G
    # (1 - 0.9534, 0.11),  # G - GHIJ
    # (1 - 0.9534, 0.06),  # G - EFGHIJ

    (0, 0.121),  # H - H
    # (1 - 0.9953, 0.12), # H - HIJ
    # (1 - 0.9966, 0.11),  # H - GHIJ
    # (1 - 0.9973, 0.06),  # H - EFGHIJ

    (0, 0.122),  # I - I
    (1 - 1, 0.12),  # I - HIJ
    (1 - 0.999, 0.11),  # I - GHIJ
    # (1 - 0.9983, 0.06), # I - EFGHIJ

    (0, 0.122),  # J - J
    # (1 - 0.9966, 0.12),  # J - HIJ
    # (1 - 0.9979, 0.11),  # J - GHIJ
    (1 - 0.9979, 0.06),  # J - EFGHIJ
]
# (((A:0.001,B:0.001)AB:0.02,
# (C:0.001,D:0.002)CD:0.01)
# ABCD:0.1,((E:0.002,F:0.002)EF:0.06,(G:0.02,(H:0.001,(I:0.001,J:0.001)IJ:0.001)HIJ:0.01)GHIJ:0.05)EFGHIJ:0.06);

x = [d[0] for d in data]
y = [d[1] for d in data]

ax.scatter(x, y, c="red")
plt.xlabel("1 - minComp")
plt.ylabel("Evolution tree distance")
plt.title("Evolution tree vs. Pangtree result - comparison")
plt.savefig("tree_distances2_tree3_nocycles500_p_1_stop_0999.png")
plt.show()
