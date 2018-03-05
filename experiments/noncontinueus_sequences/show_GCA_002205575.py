#with open('mycoplasma2_bez_nucl.maf') as inp:
with open('alignment.maf') as inp:
	a = inp.readlines()
	starts = []
	for line in a:
		line_parts = line.split()
		if len(line_parts) in [0, 1] or line_parts[1] != 'GCA_002205575':
			continue
		starts.append((int(line_parts[2]), int(line_parts[3]), int(line_parts[2]) + int(line_parts[3])))

sorted_starts = sorted(starts, key = lambda x: x[0])
for i, s in enumerate(sorted_starts):
	if i > 0 and s[0] != sorted_starts[i-1][2]:
		a = "*"
	else:
		a = ""
	print(a + str(s[0]) + " + " + str(s[1]) + " = " + str(s[2]))
		
