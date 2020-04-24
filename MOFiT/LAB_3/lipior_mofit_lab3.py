#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
import seaborn as sns


#Rysowanie wykresów
draw_chart = (len(sys.argv)>1)

#-------------------- ZMIENNE GLOBALNE --------------#
c = 1
N = 101
dx = 0.01
dt = 0.005
x = numpy.arange(0, dx*N, dx)

### Warunki początkowe
u_0 = numpy.exp(-100 * (x - 0.5)**2)
v_0 = x*0

#-------------------- KLASY	-------------------------#


#-------------------- ZADANIE 1 ---------------------#

class struna:

	def __init__(self, t_s, t_k, u_0 = u_0, v_0 = v_0, dx = dx, dt = dt, N = N):
		self.u_x_t = []
		self.v_x_t = []
		self.N = N
		self.dt = dt
		self.dx = dx
		self.u_0 = u_0
		self.v_0 = v_0
		self.t_s = t_s
		self.t_k = t_k
		self.t = numpy.arange(0, self.t_k + self.dt, self.dt)
		self.u_x_t.append(self.u_0)
		self.v_x_t.append(self.v_0)

	def cal_u_x_nt(self, u_x_t, v_x_t, a_x_t, dt):
		return u_x_t + dt*v_x_t + 0.5*a_x_t*dt**2

	def cal_v_x_nt(self, v_x_t, a_x_t, a_x_nt, dt):
		return v_x_t + 0.5*dt*(a_x_t + a_x_nt)

	def cal_a_x_t(self, u_x_t, u_nx_t, u_bx_t, dx):
		return (u_nx_t + u_bx_t - 2*u_x_t)/(dx**2)

	def wypisz(self):
		return (self.t, self.u_x_t, self.v_x_t)

	def rysuj(self, nazwa, v_min = -1, v_max = 1, plot = True, save = True):
		
		x_ticks = numpy.around(numpy.arange(self.t_s, self.t_k + 1, 1), 6)
		y_ticks = numpy.around(numpy.arange(0,1.1, 0.1), 6)

		ax = sns.heatmap(numpy.rot90(self.u_x_t), xticklabels = x_ticks, yticklabels = y_ticks[::-1], cmap = "seismic", vmin = v_min, vmax = v_max)

		ax.set_xticks(x_ticks*1/self.dt)
		ax.set_yticks(y_ticks*1/self.dx)
		ax.grid()
		ax.collections[0].colorbar.set_label("u")
		ax.figure.axes[-1].yaxis.label.set_size(18)

		ax.set(ylabel = 'x', xlabel = 't', title = 'u(x,t)')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)
		
		if plot:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + ".png")

	def energia_srednia(self,ts,tk):

		
		def helper_u(l):

			du = [0]
			for i in range(1, len(l) - 1):
				d = (l[i+1] - l[i-1])/(self.dx*2)
				du.append(d)
			du.append(0)
			return du

		E_sr = 0

		s_i = int(ts/tk * len(self.t))
		k_i = len(self.t)
		for i in range(s_i, k_i):
			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]
			E_V = 0
			E_U = 0
			du = helper_u(u_x)

			for k in range(self.N):
				E_v = v_x[k]**2 * self.dx
				E_V += E_v
				E_u = du[k]**2 * self.dx
				E_U += E_u

			E = 0.5*E_V + 0.5*E_U
			E_sr += E

		return E_sr/(self.t_k - self.t_s)


	def verleta_sztywne(self, val = 0):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero - fragment struny nie wykonuje ruchu
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		self.u_x_t[0][0] = self.u_x_t[0][-1] = val
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_x = numpy.array(helper_a(u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_x, self.dt))

			a_x_nt = numpy.array(helper_a(u_x_nt))

			v_x_nt = numpy.array(self.cal_v_x_nt(v_x, a_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)


	def verleta_luzne(self):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero, wyłącznie w celu zachowania długości list
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		self.u_x_t[0][0] = self.u_x_t[0][1]
		self.u_x_t[0][-1] = self.u_x_t[0][-2]
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_x = numpy.array(helper_a(u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_x, self.dt))

			u_x_nt[0] = u_x_nt[1]
			u_x_nt[-1] = u_x_nt[-2]

			a_x_nt = numpy.array(helper_a(u_x_nt))

			v_x_nt = numpy.array(self.cal_v_x_nt(v_x, a_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)		
		
		

#-------------------- ZADANIE 2 ---------------------#

class struna_tlumiona(struna):

	def __init__(self, t_s, t_k, beta, u_0 = u_0, v_0 = v_0, dx = dx, dt = dt, N = N):
		super().__init__(t_s, t_k, u_0 , v_0 , dx , dt , N)
		self.beta = beta

	def cal_v_x_nt(self, v_x_t, a_beta_x_t, a_x_nt, dt):
		return (v_x_t + 0.5*dt*(a_beta_x_t + a_x_nt))/(1 + self.beta * self.dt)

	def cal_a_beta_x_t(self, v_x_t, u_x_t, u_nx_t, u_bx_t, dx):
		return (u_nx_t + u_bx_t - 2*u_x_t)/(dx**2) - 2*self.beta*v_x_t

	def verleta_sztywne(self, val = 0):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero - fragment struny nie wykonuje ruchu
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		def helper_a_beta(v_x, l):
			#Przyspieszenie na brzegach przyjmujemy jako zero - fragment struny nie wykonuje ruchu
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_beta_x_t(v_x[k], l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x


		self.u_x_t[0][0] = self.u_x_t[0][-1] = val
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_beta_x = numpy.array(helper_a_beta(v_x, u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_beta_x, self.dt))

			a_x_nt = numpy.array(helper_a(u_x_nt))

			v_x_nt = numpy.array(self.cal_v_x_nt(v_x, a_beta_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)


	def verleta_luzne(self):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero, wyłącznie w celu zachowania długości list
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		def helper_a_beta(v_x, l):
			#Przyspieszenie na brzegach przyjmujemy jako zero, wyłącznie w celu zachowania długości list
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_beta_x_t(v_x[k], l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		self.u_x_t[0][0] = self.u_x_t[0][1]
		self.u_x_t[0][-1] = self.u_x_t[0][-2]
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_beta_x = numpy.array(helper_a_beta(v_x, u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_beta_x, self.dt))

			u_x_nt[0] = u_x_nt[1]
			u_x_nt[-1] = u_x_nt[-2]

			a_x_nt = numpy.array(helper_a(u_x_nt))

			v_x_nt = numpy.array(self.cal_v_x_nt(v_x, a_beta_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)		
	

#-------------------- ZADANIE 3 ---------------------#

class struna_wymuszona(struna_tlumiona):

	def __init__(self, t_s, t_k, beta, omega, x0,  u_0 = u_0, v_0 = v_0, dx = dx, dt = dt, N = N):
		super().__init__(t_s, t_k, beta, u_0 , v_0 , dx , dt , N)
		self.omega = omega
		self.x0 = x0
		self.i = int(x0 * self.N)
		self.a_F_0 = numpy.zeros(self.N)



	def cal_a_F_t(self, t):
		a = numpy.zeros(self.N)
		a[self.i] = numpy.cos(t*self.omega)
		return a

	def cal_v_x_nt(self, a_F_t, a_F_nt, v_x_t, a_beta_x_t, a_x_nt, dt):
		return (v_x_t + 0.5*dt*(a_beta_x_t + a_x_nt + a_F_t + a_F_nt))/(1 + self.beta * self.dt)

	def verleta_sztywne(self, val = 0):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero - fragment struny nie wykonuje ruchu
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		def helper_a_beta(v_x, l):
			#Przyspieszenie na brzegach przyjmujemy jako zero - fragment struny nie wykonuje ruchu
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_beta_x_t(v_x[k], l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x


		self.u_x_t[0][0] = self.u_x_t[0][-1] = val
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_F_t = self.cal_a_F_t(self.t[i])
			a_beta_x = numpy.array(helper_a_beta(v_x, u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_beta_x + a_F_t, self.dt))

			a_x_nt = numpy.array(helper_a(u_x_nt))


			a_F_nt = self.cal_a_F_t(self.t[i] + self.dt)

			v_x_nt = numpy.array(self.cal_v_x_nt(a_F_t, a_F_nt, v_x, a_beta_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)

	def verleta_luzne(self):

		def helper_a(l):
			#Przyspieszenie na brzegach przyjmujemy jako zero, wyłącznie w celu zachowania długości list
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_x_t(l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x

		def helper_a_beta(v_x, l):
			#Przyspieszenie na brzegach przyjmujemy jako zero, wyłącznie w celu zachowania długości list
			a_x = [0]
			for k in range(1, len(l) - 1):
				a_x.append(self.cal_a_beta_x_t(v_x[k], l[k], l[k+1], l[k-1], self.dx))
			a_x.append(0)
			return a_x


		self.u_x_t[0][0] = self.u_x_t[0][1]
		self.u_x_t[0][-1] = self.u_x_t[0][-2]
		
		
		for i in range(len(self.t)):

			u_x = self.u_x_t[i]
			v_x = self.v_x_t[i]

			a_F_t = self.cal_a_F_t(self.t[i])
			a_beta_x = numpy.array(helper_a_beta(v_x, u_x))
			u_x_nt = numpy.array(self.cal_u_x_nt(u_x, v_x, a_beta_x + a_F_t, self.dt))

			u_x_nt[0] = u_x_nt[1]
			u_x_nt[-1] = u_x_nt[-2]

			a_x_nt = numpy.array(helper_a(u_x_nt))


			a_F_nt = self.cal_a_F_t(self.t[i] + self.dt)

			v_x_nt = numpy.array(self.cal_v_x_nt(a_F_t, a_F_nt, v_x, a_beta_x, a_x_nt, self.dt))

			self.u_x_t.append(u_x_nt)
			self.v_x_t.append(v_x_nt)


#-------------------- ZADANIE 4 ---------------------#

class energia_srednia_struny:

	def __init__(self, t_s, t_k, beta, omega_s, omega_k, N_omega, x0, u_x_0, dt = dt, dx = dx):
		self.t_s = t_s
		self.t_k = t_k
		self.dt = dt
		self.beta = beta
		self.omega_s = omega_s
		self.omega_k = omega_k
		self.N_omega = N_omega
		self.x0 = x0
		self.dx = dx
		self.omega = numpy.linspace(omega_s, omega_k, N_omega)
		self.struny = [struna_wymuszona(t_s, t_k, beta, i, x0, u_x_0) for i in self.omega]
		self.energia_srednia = []

	def cal_energia_srednia(self):
		for i in self.struny:
			i.verleta_sztywne()
			self.energia_srednia.append(i.energia_srednia(self.t_s, self.t_k))
	
	def wypisz(self):
		return self.energia_srednia

	def rysuj(self, nazwa, plot = True, save = True):
		
		fig, ax = plt.subplots()

		x_ticks = numpy.around(numpy.arange(self.omega_s/numpy.pi, (1.01*self.omega_k)/numpy.pi , 1), 0)
		ax.plot(self.omega, self.energia_srednia)

		ax.set_xticklabels(x_ticks)
		ax.set_xticks(x_ticks * numpy.pi)

		ax.grid()

		ax.set(xlabel = '$\omega \cdot \pi$ ', ylabel = '<E>', title = '$<E>(\omega)$')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)
		
		if plot:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + ".png")
		

#-------------------- ZADANIE 1 ---------------------#

def zad_1():
	sztywne = struna(0,5)
	sztywne.verleta_sztywne()
	sztywne.rysuj('zad_1_1')

	luzne = struna(0,5)
	luzne.verleta_luzne()
	luzne.rysuj('zad_1_2')



#-------------------- ZADANIE 2 ---------------------#


def zad_2():
	a = struna_tlumiona(0, 5, 0.5)
	a.verleta_sztywne()
	a.rysuj('zad_2_1')

	b = struna_tlumiona(0, 5, 2)
	b.verleta_sztywne()
	b.rysuj('zad_2_2')

	c = struna_tlumiona(0, 5, 4)
	c.verleta_sztywne()
	c.rysuj('zad_2_3')

#-------------------- ZADANIE 3 ---------------------#

def zad_3():
	u_x_0 = numpy.zeros(N)

	a = struna_wymuszona(0, 10, 1, numpy.pi/2, 0.5, u_x_0)
	a.verleta_sztywne()
	a_l = a.wypisz()[1]
	a.rysuj('zad_3_1', numpy.amin(a_l), numpy.amax(a_l))


	b = struna_wymuszona(0, 50, 1, numpy.pi/2, 0.5, u_x_0)
	b.verleta_sztywne()
	b_l = b.wypisz()[1]
	b.rysuj('zad_3_2', numpy.amin(b_l), numpy.amax(b_l))
	

#-------------------- ZADANIE 4 ---------------------#

def zad_4_1():
	u_x_0 = numpy.zeros(N)
	a = energia_srednia_struny(16,20,1,0*numpy.pi,10*numpy.pi,100, 0.5, u_x_0)
	a.cal_energia_srednia()
	a.rysuj('zad_4_1')

def zad_4_2():
	u_x_0 = numpy.zeros(N)
	a = energia_srednia_struny(16,20,1,0,10*numpy.pi,100, 0.4, u_x_0)
	a.cal_energia_srednia()
	a.rysuj('zad_4_2')

def zad_4_3():
	u_x_0 = numpy.zeros(N)
	a = energia_srednia_struny(16,20,1,0,10*numpy.pi,100, 0.25, u_x_0)
	a.cal_energia_srednia()
	a.rysuj('zad_4_3')
#-------------------- WYKONANIE ---------------------#

#zad_1()

#zad_2()

#zad_3()

zad_4_1()

zad_4_2()

zad_4_3()



