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


n = 30
m = 150
j_1 = 60
j_2 = 90
drho = 0.1
dz = 0.1
IMAX = 5000
V_0 = 10
#-------------------- KLASY	-------------------------#

class siatka_potencjalu:

	def __init__(self, V_0 = V_0, n = n, m = m, j_1 = j_1, j_2 = j_2, iteracje = IMAX, drho = drho, dz = dz):
		self.__V_0 = V_0
		self.__iteracje = iteracje
		self.__lista_iteracji = numpy.linspace(0, self.__iteracje, self.__iteracje + 1)
		self.__siatka = (n, m)
		self.__j_1 = j_1
		self.__j_2 = j_2
		self.__drho = drho
		self.__dz = dz
		self.__rho = numpy.linspace(0, n*drho, n)
		self.__z = numpy.linspace(0, m*dz, m)
		self.__potencjal = numpy.empty(self.__siatka)
		self.__warunki_brzegowe()


	def __warunki_brzegowe(self):
		m = self.__siatka[1]
		n = self.__siatka[0]

		for j in range(self.__j_1):
			self.__potencjal[-1][j] = self.__V_0


		for j in range(self.__j_1 , self.__j_2):
			self.__potencjal[-1][j] = 0

		for j in range(self.__j_2 , m):
			self.__potencjal[-1][j] = self.__V_0

		for i in range(0, n-1):
			self.__potencjal[i][-1] = self.__potencjal[i][-2]

		for i in range(0, n-1):
			self.__potencjal[i][0] = self.__potencjal[i][1]

		for j in range(1, m-1):
			self.__potencjal[0][j] = self.__potencjal[1][j]

	def __relaksacja(self):
		for i in range(1, self.__siatka[0]-1):
			for j in range(1, self.__siatka[1]-1):
				p1 = self.__potencjal[i+1][j]
				p2 = self.__potencjal[i-1][j]
				p3 = self.__potencjal[i][j+1]
				p4 = self.__potencjal[i][j-1]
				a = 1/((2/(self.__drho**2))+(2/(self.__dz**2)))
				tmp_1 = (p1 + p2)/(self.__drho**2)
				tmp_2 =	1 / (i *self.__drho) * (p1 - p2)/(2*self.__drho)
				tmp_3 = (p3 + p4)/(self.__dz**2) 
				self.__potencjal[i][j] = a * (tmp_1 + tmp_2 + tmp_3)

	def oblicz(self, warunek = False, odwroc = False, calka_dzialania = True):
		for k in range(self.__iteracje):
			self.__relaksacja()
			self.__warunki_brzegowe()

			if not k %(self.__iteracje/10):
				print('Done', k)


	def rysuj(self, nazwa = '', draw = draw_chart, save = True):

		fig_0, ax_0 = plt.subplots()
		fit_0 = numpy.polyfit(self.__z[60:90], self.__potencjal[0][60:90], 2)
		y_0 = numpy.polyval(fit_0, self.__z[40:110])
		ax_0.plot(self.__z, self.__potencjal[0], color = 'purple')
		ax_0.plot(self.__z[40:110], y_0, color = 'lime', linestyle = 'dashed')
		ax_0.set(xlabel = '$z$', ylabel = '$V$')
		d_0 = mlines.Line2D([], [], label = 'Krzywa numeryczna: $V(0,z)$', color = 'purple', marker = None, linestyle = 'solid')
		d_0_f = mlines.Line2D([], [], label = 'Krzywa dopasowana: $V(0,z)$', color = 'lime', marker = None, linestyle = 'dashed')
		legend_0 = [d_0, d_0_f]
		ax_0.legend(handles = legend_0, loc = 'lower left', fontsize = 'x-large')

		ax_0.grid()
		ax_0.xaxis.label.set_size(18)
		ax_0.yaxis.label.set_size(18)
		fig_0 = plt.gcf()
		fig_0.set_size_inches(18.5, 10.5, forward = True)

		fig_z, ax_z = plt.subplots()
		zp = int((self.__j_1 + self.__j_2)/2)
		zz = numpy.rot90(self.__potencjal)[zp]
		fit_z = numpy.polyfit(self.__rho, zz, 2)
		rho_z = numpy.polyval(fit_z, self.__rho)
		ax_z.plot(self.__rho, zz, color = 'purple')
		ax_z.plot(self.__rho, rho_z, color = 'lime', linestyle = 'dashed')
		ax_z.set(xlabel = '$\\rho$', ylabel = '$V$')
		d_z = mlines.Line2D([], [], label = 'Krzywa numeryczna: $V(\\rho,z_p)$', color = 'purple', marker = None, linestyle = 'solid')
		d_z_f = mlines.Line2D([], [], label = 'Krzywa dopasowana: $V(\\rho,z_p)$', color = 'lime', marker = None, linestyle = 'dashed')
		legend_0 = [d_z, d_z_f]
		ax_z.legend(handles = legend_0, loc = 'lower left', fontsize = 'x-large')

		ax_z.grid()
		ax_z.xaxis.label.set_size(18)
		ax_z.yaxis.label.set_size(18)
		fig_z = plt.gcf()
		fig_z.set_size_inches(18.5, 10.5, forward = True)

		if draw:
			plt.show()

		if save:
			ax_0.get_figure().savefig(nazwa + "_r0.png", dpi = 200)
			ax_z.get_figure().savefig(nazwa + "_z.png", dpi = 200)



	def rysuj_potencjal(self, nazwa, draw = draw_chart, save = True):

		RHO, Z = numpy.meshgrid(self.__rho, self.__z)

		fig, ax = plt.subplots()
		cp = ax.contourf(RHO, Z, numpy.rot90(self.__potencjal)[::-1], cmap = 'plasma')
		cb = fig.colorbar(cp)
		cb.set_label('Potencjał $V$', size = 18)
		ax.grid()
		ax.set(ylabel = 'z', xlabel = '$\\rho$')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)
		
		fig_3d = plt.figure()
		ax_3d = plt.axes(projection='3d')
		sp = ax_3d.plot_surface(RHO, Z, numpy.rot90(self.__potencjal)[::-1], cmap = 'plasma')
		cbaxes = fig_3d.add_axes([0.05, 0.1, 0.03, 0.8])
		cb_3d = fig_3d.colorbar(cp, cax = cbaxes)
		cb_3d.set_label('Potencjał $V$', size = 18)
		ax_3d.set(ylabel = 'z', xlabel = '$\\rho$', zlabel = 'Potencjał $V$')
		ax_3d.xaxis.label.set_size(18)
		ax_3d.yaxis.label.set_size(18)
		ax_3d.zaxis.label.set_size(18)
		ax_3d.title.set_size(18)
		ax_3d.get_figure().set_size_inches(18.5, 18.5, forward = False)


		if draw:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + "_contour.png", dpi = 200)
			ax_3d.get_figure().savefig(nazwa + "_3d.png", dpi = 200)



#-------------------- ZADANIE -----------------#

def zad():
	A = siatka_potencjalu()
	A.oblicz()
	A.rysuj_potencjal('zad')
	A.rysuj()

zad()


