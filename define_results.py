from collections import defaultdict

files = open("saida.txt", "r").read().split("\n")

results = defaultdict(int)
tempos = defaultdict(int)
qtds = defaultdict(int)
print("tempo result")
for f in files[:-1]:
	V, E, L, k = map(int, f.split("/")[1].split("_")[:-1])
	str_saida = (V, E, L, k)
	results[str_saida] += int(f.split(":")[-1])
	tempos[str_saida] += int(f.split(":")[-2])
	qtds[str_saida] += 1
for k in sorted(results.keys(), key = lambda x: (x[0], x[2], x[3])):
	print(k, ":", tempos[k]/qtds[k], results[k]/qtds[k], qtds[k])
