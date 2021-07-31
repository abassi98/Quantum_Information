import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp
import math

N=input("Insert N \n")
N_lam=input("Insert N_lam \n")
k=input("How many levels?\n")
subprocess.call("./debugging.out",shell=True)

start=time.time()

var=input("Do you  want to perform RSRG? \n")
if var=='y':
	out=open("temp/parameters.txt",'w')
	out.write(N+"\n")
	out.write(N_lam+"\n")
	out.write(k+"\n")
	N_iter=input("Insert N_iter \n")
	out.write(N_iter+"\n")
	out.close()
	subprocess.call("./RSRG.out",shell=True)
var=input("Do you want to perform infinite DMRG?\n")
if var=='y':
	out=open("temp/parameters.txt",'w')
	out.write(N+"\n")
	out.write(N_lam+"\n")
	out.write(k+"\n")
	N_iter=input("Insert N_iter \n")
	out.write(N_iter+"\n")
	out.close()
	subprocess.call("./IDMRG.out",shell=True)
	
subprocess.call("gnuplot 'plot.plt'",shell=True)

finish=time.time()

print("Total computation time (s): "+ str(finish-start))
