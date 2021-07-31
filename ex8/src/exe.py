import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time
import multiprocessing as mp
import math

N=input("Insert number of sybsystems N \n")
D=input("Insert their dimension D\n")

out=open("temp/parameters.txt",'w')
out.write(N+"\n")
out.write(D+"\n")
out.close()

subprocess.call("./debugging.out",shell=True)

start=time.time()
subprocess.call("./main.out",shell=True)
finish=time.time()


print("Total computation time (s): "+ str(finish-start))
