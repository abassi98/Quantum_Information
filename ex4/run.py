import subprocess
import time

b=input("Insert N_min \n")
a=input("Insert N_max \n")
c=input("Insert N_step\n")
N_max=int(a)
N_min=int(b)
N_step=int(c)
start=time.time()
for ii in range(N_min, N_max+N_step,N_step):
	outfile=open("dimensions.txt", "w")
	outfile.write(str(ii))
	outfile.close()
	subprocess.call("./multiplication.out")
finish=time.time()
print("Total compuation time is: "+str(finish-start))
