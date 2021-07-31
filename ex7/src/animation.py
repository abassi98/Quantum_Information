import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

FFMpegWriter = animation.writers['ffmpeg']
writer = FFMpegWriter(fps=100)

matrix=np.loadtxt("data/time_ev.txt")
matrix=np.transpose(matrix)
x=matrix[0] 

#read parameters
#set default values
vec=np.loadtxt("temp/parameters.txt")
N=int(vec[0])
xmax=vec[1]
speed=vec[2]
N_time=int(vec[3])


dt=1/(speed*N_time)
fig=plt.figure()
axis=plt.axes(xlim=(-xmax/3,xmax/3),ylim=(0,1))
plt.xlabel("x")
plt.ylabel("$\psi(x,t)$")
line,=axis.plot([],[],lw=1,color='red')
line2,=axis.plot([],[],lw=1,color='blue')


def init():
	line.set_data([],[])
	line2.set_data([],[])
	return line,line2,
def V(x,i):
	t=dt*i
	return (x-speed*t)**2/2
	

def update(i):
	line.set_data(x,matrix[i+1])
	line2.set_data(x,V(x,i))
	return line,line2,
	

myanimation=animation.FuncAnimation(fig,update,init_func=init,frames=N_time,interval=5,blit=True,repeat=True)
#myanimation.save('animation.mp4',writer=writer)
plt.show()
