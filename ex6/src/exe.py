import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp
import math

N=input("Insert number of points \n")
xmax=input("Insert xmax\n")
kk=input("Insert how many eigenvectors to plot?\n")
omega=input("Insert omega\n")

outN=open("temp/N.txt", "w")
outN.write(N)
outN.close()

outxmax=open("temp/xmax.txt","w")
outxmax.write(xmax)
outxmax.close()


out=open("temp/omega.txt","w")
out.write(omega)
out.close()

out=open("temp/k.txt","w")
out.write(kk)
out.close()


subprocess.call("./debugging.out",shell=True)

start=time.time()
subprocess.call("./main.out",shell=True)
finish=time.time()

var=input("Do you want to plot?[y]/[n]")
if var=='y':
	subprocess.call("gnuplot 'plot.plt'",shell=True)

print("Total computation time: "+ str(finish-start))
