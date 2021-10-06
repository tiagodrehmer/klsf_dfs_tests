from collections import defaultdict

files = open("saida.txt", "r").read().split("\n")

results = defaultdict(int)
tempos = defaultdict(int)

for f in files:
	V, E, L, k, idd = f.split("/")[1].split("_")
	str_saida = V + "_" + E + "_" + L + "_" + k
	results[str_saida] += int(f.split(":")[-1])
	tempos[str_saida] += int(f.split(":")[-2])
for k in results.keys():
	print(k + ":", results[k], tempos[k])
