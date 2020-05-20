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

dx = 1
dy = 1

#-------------------- KLASY	-------------------------#

class siatka_potencjalu:

	def __init__(self, omega, iteracje = 10**3, nazwa = '', legend = '', kolor = 'blue', line = 'solid', pointer = None, dx = dx, dy = dy):
		self.__omega = omega
		self.__iteracje = iteracje
		self.__lista_iteracji = numpy.linspace(0, self.__iteracje, self.__iteracje + 1)
		self.__siatka = (61, 61)
		self.__dx = dx
		self.__dy = dy
		self.__x = numpy.linspace(-int(self.__siatka[0]/2), int(self.__siatka[0]/2), self.__siatka[0]//self.__dx)
		self.__y = numpy.linspace(-int(self.__siatka[1]/2), int(self.__siatka[1]/2), self.__siatka[1]//self.__dy)
		self.__nazwa = nazwa
		self.__legend = legend +'$\omega$: ' + '{:.2f}'.format(self.__omega)
		self.__kolor = kolor
		self.__line = line
		self.__pointer = pointer
		self.__gestosc_ladunku = numpy.zeros(self.__siatka)				#Ustawienie rozmiaru siatki gęstości ładunku	
		self.__ustaw_gestosc_ladunku()
		
		self.__gestosc_ladunku_z_rownania = numpy.zeros(self.__siatka) 	#Ustawienie rozmiaru siatki gęstości ładunku z odwrócenia r. Poissona
		self.__potencjal = numpy.zeros(self.__siatka)					#Poczatkowo siatka potencjalu wypelniona zerami, brzeg siatki uziemiony

		self.__S_dzialania = []

	def __ustaw_gestosc_ladunku(self):
		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):
				if i>=20 and i <=40 and j>=20 and j <=40:
					self.__gestosc_ladunku[i][j] = 1

	def __calka_dzialania(self):

		def h_pochodna_x(i, j):
			return (self.__potencjal[i+1][j] - self.__potencjal[i-1][j])/(2*self.__dx)

		def h_pochodna_y(i, j):
			return (self.__potencjal[i][j+1] - self.__potencjal[i][j-1])/(2*self.__dy)

		s = 0
		for i in range(1, self.__siatka[0]-1):
			for j in range(1, self.__siatka[1]-1):
				s+= 0.5 * (h_pochodna_x(i, j)**2 + h_pochodna_y(i, j)**2) - self.__gestosc_ladunku[i][j]*self.__potencjal[i][j]
		return s

	def __relaksacja_punktowa(self):
		for i in range(1, self.__siatka[0]-1):
			for j in range(1, self.__siatka[1]-1):
				p1 = self.__potencjal[i+1][j]
				p2 = self.__potencjal[i-1][j]
				p3 = self.__potencjal[i][j+1]
				p4 = self.__potencjal[i][j-1]
				a = (1 - self.__omega)*self.__potencjal[i][j]
				self.__potencjal[i][j] = a + (self.__omega* (p1 + p2 + p3 + p4 + (self.__gestosc_ladunku[i][j] * self.__dx * self.__dy))/4)

	def __relaksacja_globalna(self):
		tmp = numpy.zeros(self.__siatka)
		for i in range(1, self.__siatka[0]-1):
			for j in range(1, self.__siatka[1]-1):
				p1 = self.__potencjal[i+1][j]
				p2 = self.__potencjal[i-1][j]
				p3 = self.__potencjal[i][j+1]
				p4 = self.__potencjal[i][j-1]
				a = (1 - self.__omega)*self.__potencjal[i][j]
				tmp[i][j] = a + self.__omega/4 * (p1 + p2 + p3 + p4 + (self.__gestosc_ladunku[i][j] * self.__dx * self.__dy))

		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):
				self.__potencjal[i][j] = tmp[i][j]

	def __odwroc_rownanie(self):
		for i in range(1, self.__siatka[0]-1):
			for j in range(1, self.__siatka[1]-1):
				p1 = self.__potencjal[i+1][j]
				p2 = self.__potencjal[i-1][j]
				p3 = self.__potencjal[i][j+1]
				p4 = self.__potencjal[i][j-1]
				self.__gestosc_ladunku_z_rownania[i][j] = -(p1 + p2 + p3 + p4 - 4*self.__potencjal[i][j])/(self.__dx * self.__dy)

	def oblicz(self, warunek = False, odwroc = False, calka_dzialania = True):
		if not warunek:
			fun = self.__relaksacja_punktowa()
		else:
			fun = self.__relaksacja_globalna()

		for i in range(self.__iteracje + 1):
			if not warunek:
				self.__relaksacja_punktowa()
			else:
				self.__relaksacja_globalna()
			if calka_dzialania:
				self.__S_dzialania.append(self.__calka_dzialania())

			if not i%(self.__iteracje/10):
				print('Done:', i)

		if odwroc:
			self.__odwroc_rownanie()

	def wypisz(self):
		return (self.__lista_iteracji, self.__S_dzialania, (self.__kolor, self.__legend))
		
	def rysuj_gestosc_ladunku(self, nazwa, warunek = False, draw = draw_chart, save = True):
		if not warunek:
			Z = self.__gestosc_ladunku
		else:
			Z = self.__gestosc_ladunku_z_rownania

		X, Y = numpy.meshgrid(self.__x, self.__y)
		fig, ax = plt.subplots()
		cp = ax.contourf(X, Y, Z)
		cb = fig.colorbar(cp, ticks = [0,1])
		cb.ax.set_yticklabels([0, 1])
		cb.set_label('Gęstość ładunku', size = 18)
		ax.grid()
		ax.set(ylabel = 'Y', xlabel = 'X')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)
		

		fig_3d = plt.figure()
		ax_3d = plt.axes(projection='3d')
		sp = ax_3d.plot_surface(X, Y, Z, cmap = 'viridis')
		cbaxes = fig_3d.add_axes([0.05, 0.1, 0.03, 0.8])
		cb_3d = fig_3d.colorbar(cp, ticks = [0,1], cax = cbaxes)
		cb_3d.ax.set_yticklabels([0, 1])
		ax_3d.set_zticks([0,1])
		cb_3d.set_label('Gęstość ładunku', size = 18)
		ax_3d.set(ylabel = 'Y', xlabel = 'X', zlabel = 'Gęstość ładunku')
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

	def rysuj_potencjal(self, nazwa, draw = draw_chart, save = True):

		X, Y = numpy.meshgrid(self.__x, self.__y)
		fig, ax = plt.subplots()
		cp = ax.contourf(X, Y, self.__potencjal)
		cb = fig.colorbar(cp)
		cb.set_label('Potencjał $\phi$', size = 18)
		ax.grid()
		ax.set(ylabel = 'Y', xlabel = 'X')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 18.5, forward = False)
		
		fig_3d = plt.figure()
		ax_3d = plt.axes(projection='3d')
		sp = ax_3d.plot_surface(X, Y, self.__potencjal, cmap = 'viridis')
		cbaxes = fig_3d.add_axes([0.05, 0.1, 0.03, 0.8])
		cb_3d = fig_3d.colorbar(cp, cax = cbaxes)
		cb_3d.set_label('Potencjał $\phi$', size = 18)
		ax_3d.set(ylabel = 'Y', xlabel = 'X', zlabel = 'Potencjał $\phi$')
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

class wykresy_iteracji:

	def __init__(self, L, name = ''):
		self.__name = name
		self.__L = L
		self.__iteracje = [i.wypisz()[0] for i in self.__L]
		self.__S_dzialania = [i.wypisz()[1] for i in self.__L]
		self.__kolory = [i.wypisz()[2][0] for i in self.__L]
		self.__legendy = [i.wypisz()[2][1] for i in self.__L]

	def rysuj(self, nazwa, draw = draw_chart, save = True):

		#Figure
		fig, ax = plt.subplots()

		#Legend
		legend = []

		for iteracje, s_dzialania, kolor, legenda in zip(self.__iteracje, self.__S_dzialania, self.__kolory, self.__legendy):
			ax.semilogx(iteracje, s_dzialania, color = kolor)
			description = mlines.Line2D([], [], label = legenda, color = kolor, marker = None, linestyle = 'solid')
			legend.append(description)

		ax.legend(handles = legend, loc = 'upper right', fontsize = 'x-large')
		ax.set(xlabel = 'Liczba iteracji', ylabel = 'Wartość całki działania')

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if draw:
			plt.show()

		if save:
			ax.get_figure().savefig(nazwa + ".png", dpi = 200)


#Zadanie 1:
def zad_1_1():
	L = [1.85, 1.9, 1.95, 1.7, 1.5, 1.3, 1]
	C = ['blue', 'purple', 'red', 'green', 'orange', 'cyan', 'magenta']
	A = [siatka_potencjalu(i,10**3, kolor = k) for i, k in zip(L, C)]
	[i.oblicz() for i in A]
	R = wykresy_iteracji(A)
	R.rysuj('zad_1_1')

#Zadanie 1;2:
def zad__12():
	A = siatka_potencjalu(1.95, 10**3)
	A.oblicz(odwroc = True)
	A.rysuj_gestosc_ladunku('zad_1_3')
	A.rysuj_gestosc_ladunku('zad_2_1', True)
	A.rysuj_potencjal('zad_1_2')

#Zadanie 3:
def zad_3_1():
	L = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1]
	C = ['blue', 'purple', 'red', 'green', 'orange', 'cyan', 'magenta']
	A = [siatka_potencjalu(i,10**4, kolor = k) for i, k in zip(L, C)]
	[i.oblicz(True) for i in A]
	R = wykresy_iteracji(A)
	R.rysuj('zad_3_1')

def zad_3_2():	
	A = siatka_potencjalu(1.95, 10**4, kolor = 'magenta', legend = 'Relaksacja punktowa: ')
	A.oblicz()
	B = siatka_potencjalu(1, 10**4, kolor = 'lime', legend = 'Relaksacja globalna: ')
	B.oblicz()
	R = wykresy_iteracji([A,B])
	R.rysuj('zad_3_2')

#Zadanie 3 extra:
def zad_3E():
	A = siatka_potencjalu(1, 10**4)
	A.oblicz(warunek = True, odwroc = True)
	A.rysuj_gestosc_ladunku('zad_3E_2', True)
	A.rysuj_potencjal('zad_3E_3')

#zad_1_1()
#zad__12()
#zad_3_1()
#zad_3_2()
#zad_3E()



