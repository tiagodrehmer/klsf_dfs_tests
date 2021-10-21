from collections import defaultdict

files_bn = open("best_knows.txt", "r").read().split("\n")


files = open("saida.txt", "r").read().split("\n")

results_bn = defaultdict(float)
tempos_bn = defaultdict(float)
strs = defaultdict(str)
for f in files_bn:

	V, E, L, k = map(int, f.split("_")[:-1])
	keyy = (V, E, L, k)
	str_saida = str(V)+"\\_"+str(E)+"\\_"+str(L)+"\\_"+str(k)
	strs[keyy] = str_saida
	results_bn[keyy] += float(f.split(":")[-1])
	tempos_bn[keyy] += float(f.split(":")[-2])


results = defaultdict(int)
tempos = defaultdict(int)
qtds = defaultdict(int)



for f in files[:-1]:
	V, E, L, k = map(int, f.split("/")[1].split("_")[:-1])
	keyy = (V, E, L, k)
	results[keyy] += int(f.split(":")[-1])
	tempos[keyy] += int(f.split(":")[-2])
	qtds[keyy] += 1

print("\\begin{table}[]\n\\begin{tabular}{lllll}")
print("file & tempo & result & ceruli & tempo_ceruli\\\\")
for k in sorted(results.keys(), key = lambda x: (x[0], x[2], x[3])):
	print(strs[k], "&", tempos[k]/qtds[k], "&", results[k]/qtds[k], "&", str(results_bn[k]), "&", str(tempos_bn[k]), "\\\\")
print("\\end{tabular}\n\\end{table}")


