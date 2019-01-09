import os
from pathlib import Path

p1 = Path('../pang_output_01_03__10_59_42/consensus')
p2 = Path('../pang_output_01_02__21_32_27/consensus')
p3 = Path('../pang_output_01_02__21_12_37/consensus')
p4 = Path('../pang_output_01_03__12_36_27/consensus')
p5_3 = Path('../pang_output_01_03__14_54_17/consensus')
p6_3 = Path('../pang_output_01_03__14_57_48/consensus')
dir1_list = sorted(os.listdir(p2))
dir2_list = sorted(os.listdir(p4))

for i, j in zip(dir1_list, dir2_list):
    if "csv" in i:
        continue
    if i!= j:
        print("Error")
        break
    with open(p2.joinpath(i)) as i_input:
        with open(p4.joinpath(j)) as j_input:
            i_c = i_input.read()
            j_c = j_input.read()
            if i_c != j_c:
                print(f"First difference: {i}, {j}")
                break

