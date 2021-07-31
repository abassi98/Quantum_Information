import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp
import math

N_max=input("Insert N_max \n")
N_max=int(N_max)
N_min=input("Insert N_min\n")
N_min=int(N_min)
N_lam=input("Insert  N_lambda \n")
D=input("Insert their dimension D\n")
k=input("How many levels ?\n")
subprocess.call("./debugging.out",shell=True)

start=time.time()
for N in range(N_min,N_max+1):
	out=open("temp/parameters.txt",'w')
	out.write(str(N)+"\n")
	out.write(D+"\n")
	out.write(N_lam+"\n")
	out.write(k+"\n")
	out.close()
	subprocess.call("./main.out",shell=True)
finish=time.time()



print("Total computation time (s): "+ str(finish-start))
