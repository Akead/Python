#!/usr/bin/env python3
import sys
import time
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
u0 = 1

#-------------------- KLASY	-------------------------#

class siatka:
	def __init__(self, iteracje = 10**3, nazwa = '', dx = dx, dy = dy, u0 = u0):
		self.__iteracje = iteracje
		self.__siatka = (201, 102)
		self.__dx = dx
		self.__dy = dy
		self.__u0 = u0
		self.__x = numpy.arange(1, self.__siatka[0] + 1, self.__dx)

		self.__y = numpy.arange(50, self.__siatka[1] + 50, self.__dy)

		self.__fun_strumienia = numpy.zeros(self.__siatka)
		self.__potencjal =  numpy.zeros(self.__siatka)


		self.__warunki_brzegowe_strumien()
		self.__warunki_brzegowe_potencjal()

	def wypisz_strumien(self,n):
		print(numpy.rot90(self.__potencjal)[-n])

	def wypisz_potencjal(self,n):
		print(numpy.rot90(self.__potencjal)[-n])

	def __warunki_brzegowe_strumien(self):
		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):

				x = self.__x[i]
				y = self.__y[j]

				#Góra
				if y == self.__y[-1]:
					self.__fun_strumienia[i][j] = self.__u0 * y

				#Lewo i Prawo
				if y< 201 and x == self.__x[0] or y<201 and x == self.__x[-1]:
					self.__fun_strumienia[i][j] = self.__u0 * y

				#Dół

				if y == 50 and 1 < x < 201:
					self.__fun_strumienia[i][j] = self.__fun_strumienia[0][0]

				if 30 < y < 70 and 95 < x < 105:
					self.__fun_strumienia[i][j] = 'NaN'

				if  x == 95 and y < 70 or x == 105 and y < 70 or 95 <= x <= 105 and y == 70:
					self.__fun_strumienia[i][j] = self.__fun_strumienia[0][0]

	def __warunki_brzegowe_potencjal(self):
		
		#Wypełniamy tak jakby przeszkody nie było
		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):
				x = self.__x[i]
				self.__potencjal[i][j] = self.__u0 * x

		#Narzucamy warunki brzegowe
		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):

				x = self.__x[i]
				y = self.__y[j]

				#Góra
				if y == self.__y[-1]:
					self.__potencjal[i][j] = self.__u0 * x

				#Lewo
				if y< 201 and x == self.__x[0]:
					self.__potencjal[i][j] = self.__u0

				#Prawo
				if y<201 and x == self.__x[-1]:
					self.__potencjal[i][j] = self.__u0 * self.__x[-1]

				#Dół

				if y == 50 and (1 < x < 95 or 105 < x):
					self.__potencjal[i][50-50] = self.__potencjal[i][51-50]

				if (x == 95 or x == 105) and 50 < y <= 70:
					self.__potencjal[95-1][j] = self.__potencjal[94-1][j]
					self.__potencjal[105-1][j] = self.__potencjal[106-1][j]

				if y == 70 and 95 <= x < 105:
					self.__potencjal[i][70-50] = self.__potencjal[i][71-50]
				
				self.__potencjal[95-1][70-50] = (self.__potencjal[94-1][70-50] + self.__potencjal[95-1][71-50])/2
				self.__potencjal[105-1][70-50] = (self.__potencjal[106-1][70-50] + self.__potencjal[105-1][71-50])/2

				if 30 < y < 70 and 95 < x < 105:
					self.__potencjal[i][j] = 'NaN'

	def __warunki_brzegowe_potencjal_iter(self):
		for i in range(self.__siatka[0]):
			for j in range(self.__siatka[1]):

				x = self.__x[i]
				y = self.__y[j]

				if y == 50 and (1 < x < 95 or 105 < x):
					self.__potencjal[i][50-50] = self.__potencjal[i][51-50]

				if (x == 95 or x == 105) and 50 < y <= 70:
					self.__potencjal[95-1][j] = self.__potencjal[94-1][j]
					self.__potencjal[105-1][j] = self.__potencjal[106-1][j]

				if y == 70 and 95 <= x < 105:
					self.__potencjal[i][70-50] = self.__potencjal[i][71-50]
				
				self.__potencjal[95-1][70-50] = (self.__potencjal[94-1][70-50] + self.__potencjal[95-1][71-50])/2
				self.__potencjal[105-1][70-50] = (self.__potencjal[106-1][70-50] + self.__potencjal[105-1][71-50])/2

	def __relaksacja_strumien(self):

		for i in range(1, len(self.__x)-1):
			for j in range(1, self.__siatka[1]-1):

				x = self.__x[i]
				y = self.__y[j]

				if not (y <= 70 and 95 <= x <= 105):

					p1 = self.__fun_strumienia[i+1][j]
					p2 = self.__fun_strumienia[i-1][j]
					p3 = self.__fun_strumienia[i][j+1]
					p4 = self.__fun_strumienia[i][j-1]

					self.__fun_strumienia[i][j] = (p1 + p2 + p3 + p4)/4

	def __relaksacja_potencjal(self):

		for i in range(1, len(self.__x)-1):
			for j in range(1, self.__siatka[1]-1):

				x = self.__x[i]
				y = self.__y[j]

				if not (y <= 70 and 95 <= x <= 105):

					p1 = self.__potencjal[i+1][j]
					p2 = self.__potencjal[i-1][j]
					p3 = self.__potencjal[i][j+1]
					p4 = self.__potencjal[i][j-1]

					self.__potencjal[i][j] = (p1 + p2 + p3 + p4)/4

		self.__warunki_brzegowe_potencjal_iter()

	def oblicz_strumien(self):
		for k in range(self.__iteracje):

			if not(k%100):
				print(k)

			self.__relaksacja_strumien()

	def oblicz_potencjal(self):
		for k in range(self.__iteracje):

			if not(k%100):
				print(k)

			self.__relaksacja_potencjal()

	def rysuj_strumien(self, nazwa):
		
		X, Y = numpy.meshgrid(self.__x, self.__y)
		fig, ax = plt.subplots()
		V = numpy.rot90(self.__fun_strumienia)[::-1]

		C = ax.contour(X, Y, V, 25, cmap = 'viridis')
		ax.clabel(C, inline=1, fontsize=10)
		ax.fill_between([95,105],50,70, color = 'gray')
		ax.grid()
		ax.set(ylabel = 'y', xlabel = 'x', title = '$\psi(x,y)$')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)

		fig_f, ax_f = plt.subplots()
		CS = ax_f.contourf(X, Y, V, 50, cmap = 'viridis')
		ax_f.fill_between([95,105],50,70, color = 'gray')
		cb_f = fig_f.colorbar(CS)
		ax_f.grid()
		ax_f.set(ylabel = 'y', xlabel = 'x', title = '$\psi(x,y)$')
		ax_f.xaxis.label.set_size(18)
		ax_f.yaxis.label.set_size(18)
		ax_f.title.set_size(18)
		ax_f.get_figure().set_size_inches(18.5, 10.5, forward = True)


		plt.show()
		ax.get_figure().savefig(nazwa + "_contour.png", dpi = 200)
		ax_f.get_figure().savefig(nazwa + "_contour_f.png", dpi = 200)

	def rysuj_potencjal(self, nazwa):
		
		X, Y = numpy.meshgrid(self.__x, self.__y)
		fig, ax = plt.subplots()
		V = numpy.rot90(self.__potencjal)[::-1]

		z = [i for i in range(0,201,5)]
		C = ax.contour(X, Y, V, levels = z, cmap = 'viridis')
		ax.clabel(C, inline=1, fontsize=10)
		ax.fill_between([95,105],50,70, color = 'gray')


		ax.grid()
		ax.set(ylabel = 'y', xlabel = 'x', title = '$\phi(x,y)$')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)

		fig_f, ax_f = plt.subplots()
		CS = ax_f.contourf(X, Y, V, levels = 50, cmap = 'viridis')
		ax_f.fill_between([95,105],50,70, color = 'gray')
		cb_f = fig_f.colorbar(CS)
		ax_f.grid()
		ax_f.set(ylabel = 'y', xlabel = 'x', title = '$\phi(x,y)$')
		ax_f.xaxis.label.set_size(18)
		ax_f.yaxis.label.set_size(18)
		ax_f.title.set_size(18)
		ax_f.get_figure().set_size_inches(18.5, 10.5, forward = True)


		plt.show()
		ax.get_figure().savefig(nazwa + "_contour.png", dpi = 200)
		ax_f.get_figure().savefig(nazwa + "_contour_f.png", dpi = 200)



def fun_1(N):
	A = siatka(N)
	t = time.time()
	A.oblicz_strumien()
	#A.oblicz_strumien(20)
	print((time.time()-t)/60)
	A.rysuj_strumien('strumien')

def fun_2(N):
	A = siatka(N)
	t = time.time()
	#A.wypisz_potencjal(20)
	A.oblicz_potencjal()
	print((time.time()-t)/60)
	A.rysuj_potencjal('potencjal')

fun_1(10**4)
#fun_2(2*10**4)







