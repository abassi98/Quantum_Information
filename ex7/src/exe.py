import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp
import math

n=input("Insert n \n")
xmax=input("Insert xmax\n")
speed=input("Insert speed=1/T\n")
N_time=input("Insert N_time\n")
out=open("temp/parameters.txt",'w')
out.write(n+"\n")
out.write(xmax+"\n")
out.write(speed+"\n")
out.write(N_time+"\n")
out.close()

subprocess.call("./debugging.out",shell=True)

start=time.time()
subprocess.call("./main.out",shell=True)
finish=time.time()

var=input("Do you want to plot?[y]/[n]")
if var=='y':
	subprocess.call("gnuplot 'plot.plt'",shell=True)
var=input("Do you want to start the animation?[y]/[n]")
if var=='y':
	subprocess.call("python3 'animation.py'",shell=True)

print("Total computation time (s): "+ str(finish-start))
