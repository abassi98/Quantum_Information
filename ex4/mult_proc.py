import multiprocessing as mp
import subprocess
import time

b=input("Insert N_min \n")
a=input("Insert N_max \n")
c=input("Insert N_step\n")
N_max=int(a)
N_min=int(b)
N_step=int(c)

def f(ii,jj):
	name="dimensions"+str(jj)+".txt"
	exe="./multiplication"+str(jj)+".out"
	out=open(name,'w')
	out.write(str(ii))
	out.close()
	subprocess.call(exe)
	
start=time.time()	
for kk in range(N_min, N_max,4*N_step):
	processes=[mp.Process(target=f,args=(kk+(jj-1)*N_step,jj,)) for jj in range(1,5)] 
	[process.start() for process in processes]
	[process.join() for process in processes]

finish=time.time()
print('Total computation time (s): ' +str(finish-start))
var=input("Do you want to plot and fit the results?[y]/[n]")
if var=='y':
	subprocess.call("gnuplot 'fit.plt'",shell=True)




