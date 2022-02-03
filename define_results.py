from collections import defaultdict
import numpy as np
import sys
#import matplotlib.pyplot as plt

files = open(sys.argv[1], "r").read().split("\n")

#bests = open("best_knows.txt", "r").read().split("\n")


#bests_result = defaultdict(int);

#for line in bests:
#	f, t, b = line.split(":")
#	V, E, L, k, _ = map(float, f.split("_"))
#	str_saida = (V, L, k)
#	bests_result[str_saida] = float(b);


results = defaultdict(list)
tempos = defaultdict(list)
qtds = defaultdict(int)
print("tempo result")
for f in files[:-1]:
	V, E, L, k, idd = map(int, f.split("/")[1].split(":")[0].split("_"))
	str_saida = (V, L, k, idd)
	str_tempo = (V, L, k)
	results[str_saida].append(int(f.split(":")[-2]))
	tempos[str_tempo].append(int(f.split(":")[-3]))
	qtds[str_saida] += 1

print(qtds[str_saida])
groups = []
means = []
devs = []
values_per_group = defaultdict(list)
devs_per_group = defaultdict(list)
nome_per_group = defaultdict(list)

bests = []
maxs = []

for i, k in enumerate(sorted(results.keys(), key = lambda x: (x[0], x[1], x[2]))):
	group = k[:-1]
	values_per_group[group] += results[k]
	devs_per_group[group].append(np.std(results[k]))
	nome_per_group[group].append(k)




for i, g in enumerate(values_per_group.keys()):
	print("============================================")
	print(g, "{:.2f}".format(np.mean(tempos[g])),'\t', "{:.2f}".format(np.mean(values_per_group[g])), '\t', "{:.2f}".format(np.mean(devs_per_group[g])))
	print("============================================")
	for i, f in enumerate(nome_per_group[g]):
		print(f[3], '\t\t\t', "{:.2f}".format(np.mean(results[f])),'\t', "{:.2f}".format(np.std(results[f])), '\t', len(results[f]), '\t', len(results[f]) - results[f].count(min(results[f])))	
	means.append(np.mean(values_per_group[g]))
	devs.append(np.mean(devs_per_group[g]))
#	bests.append(bests_result[g])
	groups.append(str(g))




# plt.scatter(groups, bests, color='red')
# plt.errorbar(groups, means, devs, linestyle='None', marker='^')
# plt.xticks(rotation=60)
# plt.show()
