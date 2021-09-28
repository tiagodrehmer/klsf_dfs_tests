import os
import shutil
from collections import defaultdict


os.remove("klsf_instances/klsf_instances/inst200/.DS_Store")	
list_remove = os.listdir("klsf_instances")


ids = defaultdict(int)

dirr = list_remove[1]

list_dir2 = os.listdir("klsf_instances" + "/" + dirr)

for dirr2 in list_dir2[:-2]:
	files = os.listdir("klsf_instances" + "/" + dirr + "/" + dirr2)
	for f in files:
		print("klsf_instances" + "/" + dirr + "/" + dirr2 + "/" + f)
		V, E, L, k = open("klsf_instances" + "/" + dirr + "/" + dirr2 + "/" + f).read().split("\n")[0].split(" ")
		to_str =  V + "_" + E + "_" + L + "_" + k + "_"
		if ids[to_str] >= 0:
			ids[to_str] += 1
		else:
			ids[to_str] = 1
		idd = ids[to_str]
		to_str +=  str(idd)
		print(to_str)
		idd = (idd + 1) % 9
		os.rename("klsf_instances" + "/" + dirr + "/" + dirr2 + "/" + f, "klsf_instances/" + to_str)


shutil.rmtree("klsf_instances/" + list_remove[0])
shutil.rmtree("klsf_instances/" + list_remove[1])