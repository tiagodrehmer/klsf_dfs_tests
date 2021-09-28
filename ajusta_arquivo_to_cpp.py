import os

to_output = ""
for dirr in sorted(os.listdir("klsf_instances")):
	to_output += "klsf_instances/" + dirr + "\n"

f = open("lista_arquivos.txt", "w")
f.write(to_output)
f.close()