#!/usr/bin/env python3
import sys
import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import random

#Rysowanie wykresów
draw_chart = (len(sys.argv)<=1)

#------------------------ KLASY	---------------------#


#----------------------- ZADANIE 1 ------------------#

class Schemat_Metropolisa:

	def __init__(self, x_0, N):
		self.__x_0 = x_0
		self.__N = N

		self.__x_1 = 0
		self.__x_2 = 0
		self.__x_3 = 0
		self.__x_4 = 0

	def __prawdopodobienstwo(self, x):
		return numpy.pi**(-0.5) * numpy.exp(-x**2)
	

	def oblicz(self):
		i = 0
		x = self.__x_0

		while(i != self.__N):
			dx = random.uniform(-0.25, 0.25)
			x_n = x + dx
			warunek = random.random()
			if warunek < self.__prawdopodobienstwo(x_n)/self.__prawdopodobienstwo(x):
				self.__x_1 += x**1
				self.__x_2 += x**2
				self.__x_3 += x**3
				self.__x_4 += x**4
				x = x_n
				i += 1
		return (self.__x_1/self.__N, self.__x_2/self.__N, self.__x_3/self.__N, self.__x_4/self.__N)

class Lista_Metropolisa:

	def __init__(self, L, expected, legend_labels, colors, name, x_0 = 0):
		self.__L = L
		self.__x_0 = x_0

		self.__XN = []
		self.__expected = expected
		self.__legend_labels = legend_labels
		self.__colors = colors
		self.__name = name

	def oblicz(self):
		for i in self.__L:
			A = Schemat_Metropolisa(self.__x_0, i)
			self.__XN.append(A.oblicz())
			print('Obliczone:\t',i)

	def print(self, draw = draw_chart, save = True):
		self.oblicz()
		#Figure
		fig, ax = plt.subplots()
		#Legend
		legend = []

		for i, color, legend_label, expect in zip(numpy.rot90(self.__XN)[::-1], self.__colors, self.__legend_labels, self.__expected):
			ax.semilogx(self.__L, i, color = color, marker = '*')
			ax.axhline(expect, color = color, linestyle = 'dashed', alpha = 0.5)
			description = mlines.Line2D([], [], label = legend_label, color = color)
			legend.append(description)

		ax.legend(handles = legend, loc = 'upper right', fontsize = 'x-large')

		#Label XY
		ax.set(ylabel = '$I_n(l)$', xlabel = 'l')

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if save:
			fig.savefig(self.__name+'.png', dpi = 100)

		if draw:
			plt.show()


#----------------------- ZADANIE 2 ------------------#

class Schemat_Metropolisa_2D:

	def __init__(self, x_0, y_0, N, warunek = False, name = '', linestyle = 'None'):
		self.__x_0 = x_0
		self.__y_0 = y_0
		self.__N = N
		self.__warunek = warunek
		self.__linestyle = linestyle
		self.__name = name
	
		self.__V = 0

		self.__X_H = []
		self.__Y_H = []

		self.__X_M = []
		self.__Y_M = []


	def __prawdopodobienstwo(self, x, y):
		return abs((numpy.pi**(-0.5) * numpy.exp(-0.5*(x**2 + y**2)))**2)
	

	def oblicz(self):
		i = 0
		x = self.__x_0
		y = self.__y_0

		while(i != self.__N):
			dx = random.uniform(-0.25, 0.25)
			dy = random.uniform(-0.25, 0.25)
			x_n = x + dx
			y_n = y + dy

			warunek = random.random()
			if warunek < self.__prawdopodobienstwo(x_n, y_n)/self.__prawdopodobienstwo(x, y):
				self.__V += 0.5 * (x**2 + y**2)
				x = x_n
				y = y_n
				if self.__warunek:
					self.__X_H.append(x_n)
					self.__Y_H.append(y_n)
				i += 1
			elif self.__warunek:
				self.__X_M.append(x_n)
				self.__Y_M.append(y_n)


		return(self.__V/self.__N)

	def plot(self, draw = draw_chart, save = True):
		if self.__warunek:
			self.oblicz()

			fig, ax = plt.subplots()

			#Missed
			ax.plot(self.__X_M, self.__Y_M, marker = '.', color = 'deepskyblue', linestyle = 'None')

			#Hitted
			ax.plot(self.__X_H, self.__Y_H, '.', alpha = 0.4, color = 'orange', linestyle = self.__linestyle)

			#Legend
			s_M = "Miss:" +'  '+ str(len(self.__X_M))
			s_H = "Hit:" +'   '+ str(len(self.__X_H))
			desc_M = mlines.Line2D([], [], label = s_M, color = 'deepskyblue', linestyle = 'None', marker = '.')
			desc_H = mlines.Line2D([], [], label = s_H, color = 'orange', linestyle = 'None', marker = '.')
			legend = [desc_M, desc_H]

			ax.legend(handles = legend, loc = 'upper right', fontsize = 'x-large')
			ax.grid()

			#Label XY
			ax.set(ylabel = 'x', xlabel = 'y')

			ax.xaxis.label.set_size(18)
			ax.yaxis.label.set_size(18)
			fig = plt.gcf()
			fig.set_size_inches(10.5, 10.5, forward = True)

			if save:
				fig.savefig(self.__name+'.png', dpi = 100)

			if draw:
				plt.show()

class Lista_Metropolisa_2D:

	def __init__(self, L, expected, name, x_0 = 0, y_0 = 0):
		self.__L = L
		self.__x_0 = x_0
		self.__y_0 = y_0

		self.__V = []
		self.__expected = expected
		self.__name = name

	def oblicz(self):
		for i in self.__L:
			A = Schemat_Metropolisa_2D(self.__x_0, self.__y_0, i)
			self.__V.append(A.oblicz())
			print('Obliczone:\t',i)

	def print(self, draw = draw_chart, save = True):
		self.oblicz()
		#Figure
		fig, ax = plt.subplots()
		
		ax.semilogx(self.__L, self.__V, color = 'orange', marker = '*', markeredgecolor = 'blue', markerfacecolor = 'blue')
		ax.axhline(self.__expected, color = 'red', linestyle = 'dashed')

		#Label XY
		ax.set(ylabel = '$<E_p> = < \\frac{x^2 + y^2}{2} >$', xlabel = 'l')

		#Legend
		s_E = "Wartość spodziewana:" +'  '+ str(self.__expected)
		s_C = "Dane doświadczalne"
		desc_E = mlines.Line2D([], [], label = s_E, color = 'red', linestyle = 'dashed')
		desc_C = mlines.Line2D([], [], label = s_C, color = 'orange', marker = '*', markeredgecolor = 'blue', markerfacecolor = 'blue')
		legend = [desc_E, desc_C]

		ax.legend(handles = legend, loc = 'upper left', fontsize = 'x-large')

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if save:
			fig.savefig(self.__name+'.png', dpi = 100)

		if draw:
			plt.show()


#-------------------- ZADANIE 1 ---------------------#

def zad_1(a, b, name):
	
	#Przygotowanie danych od 10^a do 10^b
	l = []
	for i in range(a,b):
		for k in range(10**i, 10**(i+1), 10**(i)):
			l.append(k)
	l.append(10**b)
	
	colors = ['blue', 'red', 'orange', 'lime']
	legend_labels = ['n = 1', 'n = 2', 'n = 3', 'n = 4']
	expected = [0, 0.5, 0, 0.75]

	a = Lista_Metropolisa(l, expected, legend_labels, colors, name)
	a.print()


#zad_1(0, 7, 'zad_test')

#-------------------- ZADANIE 2 ---------------------#

def zad_2_1():
	A = Schemat_Metropolisa_2D(0, 0, 5*10**1, True, 'zad_2_1_1', 'solid')
	B = Schemat_Metropolisa_2D(0, 0, 10**2, True, 'zad_2_1_2', 'solid')
	C = Schemat_Metropolisa_2D(0, 0, 5*10**2, True, 'zad_2_1_2', 'solid')
	D = Schemat_Metropolisa_2D(0, 0, 10**3, True, 'zad_2_1_3', 'solid')
	E = Schemat_Metropolisa_2D(0, 0, 10**4, True, 'zad_2_1_4')
	F = Schemat_Metropolisa_2D(0, 0, 10**5, True, 'zad_2_1_5')
	l = [A, B, C, D, E, F]
	for i in l:
		i.plot()

def zad_2_2(a, b, name):
	
	#Przygotowanie danych od 10^a do 10^b
	l = []
	for i in range(a,b):
		for k in range(10**i, 10**(i+1), 10**(i)):
			l.append(k)
	l.append(10**b)
	
	expected = 0.5

	a = Lista_Metropolisa_2D(l, expected, name)
	a.print()


#zad_2_1()

zad_2_2(0, 7, 'zad_2_2')






