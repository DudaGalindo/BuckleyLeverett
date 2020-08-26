import numpy as np
from read_inputs import data
from matplotlib import pyplot as plt
from matplotlib import animation

global Sw
global n
datas = np.load('results_Buckley_Leverett_FOUM.npy', allow_pickle=True)
Sw = np.empty((datas.shape[0]-1,data['mesh_elements']))
t = np.empty(datas.shape[0]-1)
n = datas.shape[0]-1
i=0
for data in datas[i+1:]:
    Sw[i] = data[1]
    t[i] = data[0]
    i+=1

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')

for i in range(n):
    Sw_plot = Sw[i]
    x = np.linspace(0,1,n)
    interval = 1
    artists = init()
    artists = func(Sw_plot)
    fig.canvas.draw_idle()
    fig.canvas.start_event_loop(interval)
# call the animator.  blit=True means only re-draw the parts that have changed.

plt.show()
