import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp

start1=0
finish1=0
def f():
	subprocess.call("./norm_spacings.out")
print("This program computes the normalized spacings for random matrices \n")
var=input("Do you want to make a simulation?[y]/[n]\n")
if var=='y':
	maximum=input("Insert the dimension\n")
	maximum=int(maximum)	
	subprocess.call("./debugging.out")
	start1=time.time()
	for ii in range(1,5):
		outdim=open("fortran/temp/dimension.txt", "w")	
		outdim.write(str(maximum))
		outdim.close()
		outdiv=open("fortran/temp/division.txt", "w")
		outdiv.write(str(maximum))
		outdiv.close()
		processes=[mp.Process(target=f) for jj in range(1,5)] 
		[process.start() for process in processes]
		[process.join() for process in processes]
	finish1=time.time()
start2=0
finish2=0
var=input("Do you want to compute the probability distributions?[y]/[n]")

if var=='y':
	nbins=input("Insert number of bins \n")
	xmin=input("Insert xmin \n")
	xmax=input("Insert xmax \n")
	outnbins=open("fortran/temp/nbins.txt", "w")
	outnbins.write(nbins)
	outnbins.close()
	outxmin=open("fortran/temp/xmin.txt", "w")
	outxmin.write(xmin)
	outxmin.close()
	outxmax=open("fortran/temp/xmax.txt", "w")
	outxmax.write(xmax)
	outxmax.close()
	start2=time.time()
	subprocess.call("./prob_distribution.out") 
	finish2=time.time()
print("Total compuation time is: "+str(finish1-start1+finish2-start2))
var=input("Do you want to plot and fit the results?[y]/[n]")

if var=='y':
	subprocess.call("gnuplot 'plot.plt'",shell=True)
	
