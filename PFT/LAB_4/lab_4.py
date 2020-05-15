#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
import seaborn as sns

#Chart drawing
draw_chart = (len(sys.argv)<=1)


#-------------------- GLOBAL VARIABLES --------------#

N = 50
dt = 0.02
delta = 0.1
alpha = 1
m = 1



#-------------------- GLOBAL FUNCTIONS ---------------#

	#METHOD RK4
def RK4(start_cond, function, dt):
	k_1 = function(start_cond)
	k_2 = function(start_cond + k_1*dt/2)
	k_3 = function(start_cond + k_2*dt/2)
	k_4 = function(start_cond + k_3*dt)
	
	final = start_cond + dt/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)

	return final


#------------------------- CLASSES -------------------#

class atoms_chain:

	def __init__(self, x_0, v_0, nt, dt = dt, N = N, delta = delta, aplha = alpha, m = m):
		self._x_0 = x_0
		self._v_0 = v_0
		self._nt = nt
		self._dt = dt
		self._N = N + 1
		self._alpha = alpha
		self._delta = delta
		self._m = m
		self._x_t = [x_0]
		self._u_t = []
		self._v_t = [v_0]
		self._t = numpy.linspace(0, nt*dt, nt + 1)
		self._tt = 0
		self._T = []
		self._U = []
		self._E = []

	def _helper_a(self, l):
		a_x = [0]
		for i in range(1, len(l)-1):
			a_x.append(self._alpha/self._m * (l[i-1] - 2*l[i] + l[i+1]))
		a_x.append(0)
		return a_x

	def _function(self, s_vector):

		s_0 = s_vector[0]
		s_1 = s_vector[1]

		d_0 = s_1
		d_0[0] = d_0[-1] = 0

		d_1 = self._helper_a(s_0)

		return numpy.array([d_0, d_1])

	def _kinetyczna(self):
		for i in self._v_t:
			self._T.append(self._m/2 * sum(i**2))

	
	def _potencjalna(self):
		for i in self._x_t:
			E = 0
			for k in range(1,len(i)):
				E += self._alpha/2 * (i[k-1] - i[k] + self._delta)**2
			self._U.append(E)

	def _calkowita(self):
		for i,k in zip(self._T, self._U):
			self._E.append(i+k)


	def calculate(self):
		l = [i*delta for i in range(N+1)]
		for i in range(self._nt):

			x_i_0 = self._x_t[i]
			v_i_0 = self._v_t[i]
			self._u_t.append(x_i_0 - l)

			new_conditions = RK4([x_i_0, v_i_0], self._function, self._dt)
			x_n = new_conditions[0]
			v_n = new_conditions[1]

			self._tt += dt
			self._x_t.append(x_n)
			self._v_t.append(v_n)

		self._potencjalna()
		self._kinetyczna()
		self._calkowita()

	def rysuj(self, nazwa, x = 10, plot = draw_chart, save = True):

		x_ticks = numpy.around(numpy.arange(0, self._dt*self._nt + 1 , x), 6)

		fig, ax = plt.subplots()

		ax.plot(self._t, self._T, color = 'red')
		ax.plot(self._t, self._U, color = 'blue')
		ax.plot(self._t, self._E, color = 'green')

		ax.set_xticklabels = x_ticks
		ax.set_xticks(x_ticks)

		#Legend

		description_T = mlines.Line2D([], [], label = '$E_k$', color = 'red', marker = None)
		description_U = mlines.Line2D([], [], label = '$E_p$', color = 'blue', marker = None)
		description_E = mlines.Line2D([], [], label = '$E_c$', color = 'green', marker = None)
		legend = [description_T, description_U, description_E]
		ax.legend(handles = legend, loc = 'lower right', fontsize = 'x-large')

		#Label XY
		ax.set(ylabel = 'E', xlabel = 't')

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if save:
			fig.savefig(nazwa+'.png', dpi = 100)

		if plot:
			plt.show()


	def rysuj_heat(self, nazwa, x = 10, y = 5, v_min = -0.05, v_max = 0.05, plot = draw_chart, save = True):


		x_ticks = numpy.around(numpy.arange(0, self._dt*self._nt + 1, x), 6)
		y_ticks = numpy.around(numpy.arange(0,self._N, y), 6)

		ax = sns.heatmap(numpy.rot90(self._u_t) , xticklabels = x_ticks, yticklabels = y_ticks[::-1], cmap = "seismic", vmin = v_min, vmax = v_max)

		ax.set_xticks(x_ticks*1/self._dt)
		ax.set_yticks(y_ticks)
		ax.grid()
		ax.collections[0].colorbar.set_label("u")
		ax.figure.axes[-1].yaxis.label.set_size(18)

		ax.set(ylabel = 'x', xlabel = 't', title = 'u(x,t)')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)
		ax.get_figure().subplots_adjust(bottom=0.14)
		if plot:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + ".png")


class resonanse_atoms_chain(atoms_chain):

	def __init__(self, n, F_0, x_0, v_0, dt = dt, N = N, delta = delta, aplha = alpha, m = m):
		self._n = n
		self._F_0 = F_0
		self._x_0 = x_0
		self._v_0 = v_0
		self._dt = dt
		self._N = N + 1
		self._alpha = alpha
		self._delta = delta
		self._m = m
		self._omega_n = 2* numpy.sqrt(self._alpha/self._m) * abs(numpy.sin(self._n*numpy.pi/(2*self._N)))
		self._t_max = 20*2*numpy.pi/self._omega_n
		self._nt = int(self._t_max/self._dt)
		self._x_t = [x_0]
		self._u_t = []
		self._v_t = [v_0]
		self._t = numpy.linspace(0, self._t_max, self._nt + 1)
		self._tt = 0
		self._T = []
		self._U = []
		self._E = []

	def _helper_a(self, l):
		a_x = [0]
		for i in range(1, len(l)-1):
			a_x.append(self._alpha/self._m * (l[i-1] - 2*l[i] + l[i+1]))
		a_x.append(0)
		a_x[1] = self._alpha/self._m * (l[0] - 2*l[1] + l[2]) + self._F_0/self._m * numpy.sin(self._omega_n * self._tt)
		return a_x


#----------------------- EXECUTION --------------------------

def zad_1():

	def x_start():
		def start_h(i):
			x_max = delta * (N+1)
			sigma = 3*delta
			return i*delta + delta/3 * numpy.exp(-((delta*i-x_max/2)**2)/(2*sigma**2))

		x = []
		for i in range(N+1):
			x.append(start_h(i))
		return x

	x_s = numpy.array(x_start())
	v_s = numpy.zeros(N+1)

	A = atoms_chain(x_s , v_s, 5000)
	A.calculate()
	A.rysuj_heat('zad_1_H')
	A.rysuj('zad_1_E')

def zad_2(n, nazwa_1, nazwa_2):

	def x_start():
		def start_h(i):
			return i*delta

		x = []
		for i in range(N+1):
			x.append(start_h(i))
		return x

	x_s = numpy.array(x_start())
	v_s = numpy.zeros(N+1)

	A = resonanse_atoms_chain(n, 0.01, x_s , v_s)
	A.calculate()
	A.rysuj_heat(nazwa_1, 400)
	A.rysuj(nazwa_2, 400)

zad_1()

#zad_2(0.9, 'zad_2_H_09', 'zad_2_E_09')
#zad_2(1, 'zad_2_H_1', 'zad_2_E_1')
#zad_2(1.1, 'zad_2_H_11', 'zad_2_E_11')
#zad_2(1.5, 'zad_2_H_15', 'zad_2_E_15')
#zad_2(2, 'zad_2_H_20', 'zad_2_E_20')
#zad_2(5, 'zad_2_H_50', 'zad_2_E_50')



