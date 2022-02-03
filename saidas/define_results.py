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
#	str_saida = (V, L)
#	bests_result[str_saida] = float(b);




results_per_seed = defaultdict(list)
results = defaultdict(list)
tempos = defaultdict(list)
qtds = defaultdict(int)
value_seeds = np.zeros(30)


mmm = 0

v = len(files[:-1])%200

print(v)
for i, f in enumerate(files[:-(v + 1)]):
	V, E, L, k, idd = map(int, f.split("/")[1].split(":")[0].split("_"))
	str_saida = (V, L, k, idd)
	str_tempo = (V, L, k)
	str_g_s = (V, L, k, idd)
	results[str_saida].append(int(f.split(":")[-2]))
	tempos[str_tempo].append(int(f.split(":")[-3]))
	if i % 200 == 0:
            mmm+=1
	value_seeds[mmm-1] += results[str_saida][-1]
	results_per_seed[(str_tempo, mmm-1)].append(results[str_saida][-1])
	qtds[str_saida] += 1

print(f)

if(i>200):
	best_seed = np.argmin(value_seeds[:int(np.floor(i/200))])
	value_b_seed = np.min(value_seeds[:int(np.floor(i/200))])
else:
	value_b_seed = 0
	best_seed = 0


groups = []
means = []
devs = []
values_per_group = defaultdict(list)
devs_per_group = defaultdict(list)
nome_per_group = defaultdict(list)
its_per_group = defaultdict(list)
bests = []
maxs = []

for i, k in enumerate(sorted(results.keys(), key = lambda x: (x[0], x[1], x[2]))):
	group = k[:-1]
	values_per_group[group] += results[k]
	devs_per_group[group].append(np.std(results[k]))
	nome_per_group[group].append(k)

media_geral_value = 0
media_geral_tempo = 0

total_not = 0
total_por_cento = 0
total_dev = 0
total_mean = 0
total_t = 0
total_bs = 0

total_not2 = 0
total_por_cento2 = 0
total_dev2 = 0
total_mean2 = 0
total_t2 = 0
total_bs2 = 0
print("f \t\t t \t v \t sd \t bs("+str(best_seed)+") \t %")
print("========================================")

#data = []


for i, g in enumerate(values_per_group.keys()):
	media_geral_value += np.mean(values_per_group[g])
	media_geral_tempo += np.mean(tempos[g])

	somat = 0
	minn = 0
	ttttt = 0	
	#data.append(np.zeros(30))
	for j, f in enumerate(nome_per_group[g]):
	#	print(j, ":", f[3], '\t\t\t\t', "{:.2f}".format(np.mean(results[f])),'\t', "{:.2f}".format(np.std(results[f])), '\t', len(results[f]), '\t', len(results[f]) - results[f].count(min(results[f])))
		total_not += len(results[f]) - results[f].count(min(results[f]))
		minn += min(results[f])
		ttttt += len(results[f])
		#data[-1] += results[f]
	minn = minn/10
	#data[-1] = data[-1]/10
	print(g, "\t{:.2f}".format(np.mean(tempos[g])),'\t', "{:.2f}".format(np.mean(values_per_group[g])), '\t', "{:.2f}".format(np.mean(devs_per_group[g])), '\t', "{:.2f}".format(np.mean(results_per_seed[(g,best_seed)])), "\t", "{:.2f}".format(100 - total_not*100/300), "\t", ttttt)
	if(i % 4 == 3):
        	print("========================================")

	if(i > 9):
		total_por_cento2 += 100 - total_not*100/300
		total_not2 = 0
		total_dev2 += np.mean(devs_per_group[g])
		total_t2 += np.mean(tempos[g])
		total_mean2 += np.mean(values_per_group[g])
		total_bs2 += np.mean(results_per_seed[(g,best_seed)])

	total_por_cento += 100 - total_not*100/300
	total_not = 0
	total_dev += np.mean(devs_per_group[g])
	total_t += np.mean(tempos[g])
	total_mean += np.mean(values_per_group[g])
	total_bs += np.mean(results_per_seed[(g,best_seed)])
	means.append(np.mean(values_per_group[g]))
	devs.append(np.mean(devs_per_group[g]))
	#if(minn< bests_result[g]):
	#	bests.append(minn)
	#else:
	#	bests.append(bests_result[g])
	groups.append(str(g))

print("average", "\t{:.2f}".format(total_t/20),'\t', "{:.2f}".format(total_mean/20), "\t", "{:.2f}".format(total_dev/20), '\t', "{:.2f}".format(total_bs/20), "\t", "{:.2f}".format(total_por_cento/20))
print("average", "\t{:.2f}".format(total_t2/10),'\t', "{:.2f}".format(total_mean2/10), "\t", "{:.2f}".format(total_dev2/10), '\t', "{:.2f}".format(total_bs2/10), "\t", "{:.2f}".format(total_por_cento2/10))
		
ss = [100 for i in range(20)]


#if (sys.argv[2] == "1"):
	#plt.scatter(groups, bests, color='red', facecolors='none', s = 100, marker = 'o')
	#plt.errorbar(groups, means, devs, linestyle='None', marker='.', markersize = 5)
	#plt.xticks(rotation=90)
	#plt.subplots_adjust(bottom=0.23)
	#plt.savefig('figure1.pdf')

	# fig = plt.figure(figsize =(10, 7))
 
	# # Creating plot
	# plt.boxplot(data)
	# plt.xticks(range(1, 21) , groups, rotation=90)
	#plt.scatter(groups, bests, color='red', s = ss)
	# # show plot
	# plt.show()
