#!/usr/bin/env python3
import sys
import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits import mplot3d
from matplotlib import cm

#Rysowanie wykresów
draw_chart = (len(sys.argv)<=1)

#-------------------- ZMIENNE GLOBALNE --------------#



#-------------------- KLASY	-------------------------#
class siatka_wezlow:

	def __init__(self, alpha, L, K = 41):
		self.__alpha = alpha
		self.__L = L
		self.__K = K
		self.__u = 1
		self.__E = 1
		self.__R = 1
		self.__sigma = 1
		self.__omega = 1
		self.__omega_x = self.__omega * numpy.sin(alpha)
		self.__omega_y = self.__omega * numpy.cos(alpha)
		self.__OMEGA = numpy.array([self.__omega_x, self.__omega_y, 0])
		self.__N = 201
		self.__M = 201
		self.__delta_x = 2*L/K
		self.__delta_y = 2*L/K
		self.__delta_th = numpy.pi/self.__N
		self.__delta_fi = 2*numpy.pi/self.__M

		self.__theta = numpy.linspace(0, numpy.pi, self.__N)
		self.__fi = numpy.linspace(0, 2*numpy.pi, self.__M)
		self.__x = numpy.linspace(-L, L, K)
		self.__y = numpy.linspace(-L, L, K)

		self.__Vr = []
		self.__Ex = []
		self.__Ey = []
		self.__Ez = []

		self.__Bx = []
		self.__By = []
		self.__Bz = []


	def __r_ij(self, i, j):
		x = self.__x[i]
		y = self.__y[j]
		z = 0
		return numpy.array([x, y, z])

	def __R_ij(self, i, j):
		R_x = self.__R * numpy.sin(self.__theta[i]) * numpy.cos(self.__fi[j])
		R_y = self.__R * numpy.sin(self.__theta[i]) * numpy.sin(self.__fi[j])
		R_z = self.__R * numpy.cos(self.__theta[i])
		return numpy.array([R_x, R_y, R_z])


	def __wsp(self, i, N):
		k = 1 if i == 0 or i == N else 4 if i%2 else 2
		return k


	def __dlugosc_vector(self, x):
		v = numpy.array(x)
		return numpy.sqrt(v.dot(v))

	def __iloczyn_vector(self, x, y):
		X = numpy.array(x)
		Y = numpy.array(y)
		return numpy.cross(X,Y)

	def __potencjal(self, i, j):
		a = self.__sigma * self.__R**2 / (4*numpy.pi * self.__E) * self.__delta_th * self.__delta_fi/9
		V = 0
		for k in range(self.__N):
			for l in range(self.__M):
				R_diff = self.__dlugosc_vector((self.__r_ij(i,j) - self.__R_ij(k,l)))
				V += self.__wsp(k, self.__N) * self.__wsp(l, self.__M) * numpy.sin(self.__theta[k])/R_diff
		return a*V

	def __pole_E(self, i,j):
		a = self.__sigma * self.__R**2 / (4*numpy.pi * self.__E) * self.__delta_th * self.__delta_fi/9
		E = [0,0,0]
		for k in range(self.__N):
			for l in range(self.__M):
				R_diff = self.__dlugosc_vector((self.__r_ij(i,j) - self.__R_ij(k,l)))
				RR_diff = self.__r_ij(i,j) - self.__R_ij(k,l)
				E += self.__wsp(k, self.__N) * self.__wsp(l, self.__M) * numpy.sin(self.__theta[k])*RR_diff/R_diff**3
		return a*E

	def __pole_B(self, i,j):
		a = -self.__sigma * self.__u * self.__R**2 / (4*numpy.pi) * self.__delta_th * self.__delta_fi/9
		B = [0,0,0]
		for k in range(self.__N):
			for l in range(self.__M):
				R_diff = self.__dlugosc_vector((self.__r_ij(i,j) - self.__R_ij(k,l)))
				RR_diff = self.__r_ij(i,j) - self.__R_ij(k,l)
				WR = self.__iloczyn_vector(self.__OMEGA, self.__R_ij(k,l))
				RD_WR = self.__iloczyn_vector(RR_diff, WR)
				B += self.__wsp(k, self.__N) * self.__wsp(l, self.__M) * numpy.sin(self.__theta[k]) * RD_WR/(R_diff**3)
		return a*B


	def oblicz(self):

		for i in range(self.__K):
			V = []

			for j in range(self.__K):
				V.append(self.__potencjal(i,j))
				E = (self.__pole_E(i,j))
				B = (self.__pole_B(i,j))
				self.__Ex.append(E[0])
				self.__Ey.append(E[1])
				self.__Ez.append(E[2])
				self.__Bx.append(B[0])
				self.__By.append(B[1])
				self.__Bz.append(B[2])
				print(i,j)

			self.__Vr.append(V)


	def rysuj_potencjal(self, nazwa, draw = draw_chart, save = True):
		X, Y = numpy.meshgrid(self.__x, self.__y)
		fig, ax = plt.subplots()
		x = [numpy.cos(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		y = [numpy.sin(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		ax.plot(x,y, color = 'black')
		cp = ax.contourf(X, Y, self.__Vr, cmap = 'plasma')
		cb = fig.colorbar(cp)
		cb.set_label('Potencjał $V$', size = 18, fontsize = 20)
		ax.grid()
		ax.tick_params(labelsize = 'large')
		ax.set(ylabel = 'y', xlabel = 'x')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)

		if draw:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + "_contour.png", dpi = 200)


	def rysuj_pole_E(self, nazwa, draw = draw_chart, save = True):
		X, Y = numpy.meshgrid(self.__x, self.__y)

		fig, ax = plt.subplots()
		x = [numpy.cos(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		y = [numpy.sin(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		ax.plot(x,y, color = 'red')
		cp = ax.quiver(Y, X, self.__Ex, self.__Ey)
		ax.grid()
		ax.tick_params(labelsize = 'large')
		ax.set(ylabel = 'y', xlabel = 'x')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)

		if draw:
			plt.show()


		if save:
			ax.get_figure().savefig(nazwa + "_E.png", dpi = 200)


	def rysuj_pole_B(self, nazwa, draw = draw_chart, save = True):
		X, Y = numpy.meshgrid(self.__x, self.__y)

		fig, ax = plt.subplots()
		x = [numpy.cos(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		y = [numpy.sin(i) for i in numpy.linspace(0, 2*numpy.pi, 1000)]
		ax.plot(x,y, color = 'red')
		cp = ax.quiver(Y, X, self.__Bx, self.__By)
		ax.grid()
		ax.tick_params(labelsize = 'large')
		ax.set(ylabel = 'y', xlabel = 'x')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)

		if draw:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + "_B.png", dpi = 200)


	def rysuj_Vy(self, nazwa, draw = draw_chart, save = True):

		fig, ax = plt.subplots()
		fig.subplots_adjust(bottom=0.15)
		ax.plot(self.__y, numpy.rot90(self.__Vr)[self.__K//2], color = 'red')
		ax.grid()
		ax.tick_params(labelsize = 'large')
		ax.set(ylabel = 'V(0,y,0)', xlabel = 'y')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)

		if draw:
			plt.show()


		if save:
			ax.get_figure().savefig(nazwa + "_y.png", dpi = 200)


def zad_1():
	A = siatka_wezlow(0,3)
	A.oblicz()
	A.rysuj_potencjal('1')
	A.rysuj_pole_E('1')
	A.rysuj_pole_B('1')
	A.rysuj_Vy('1')

def zad_2():
	A = siatka_wezlow(numpy.pi/4,3)
	A.oblicz()
	A.rysuj_potencjal('2')
	A.rysuj_pole_E('2')
	A.rysuj_pole_B('2')
	A.rysuj_Vy('2')

zad_1()
#zad_2()


