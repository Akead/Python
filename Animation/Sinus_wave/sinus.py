#!/usr/bin/env python3

import numpy
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

import matplotlib
matplotlib.use("Agg")

lx = []
ly = []
for i in range(1000):
	x = numpy.linspace(0,4,1000)
	y = numpy.sin(2 * numpy.pi * (x - 0.01 * i))
	lx.append(x)
	ly.append(y)

def init():
	line.set_data([], [])
	return line,

def animate(i):
	x = lx[i]
	y = ly[i]
	line.set_data(x, y)
	return line,

fig = plt.figure()
ax = plt.axes(xlim = (0,4), ylim = (-4,4))
line, = ax.plot([], [])

anim = FuncAnimation(fig, animate, init_func = init, frames = len(lx), interval = 20, blit = True, save_count = 1000)

plt.show()

anim.save('Sinus_wave.mp4', progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
