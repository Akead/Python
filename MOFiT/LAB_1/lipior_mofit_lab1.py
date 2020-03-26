#!/usr/bin/env python3
import sys
import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


file = open('wynik.txt','w')

#Rysowanie wykresów
draw_chart = (len(sys.argv)>1)

#-------------------- ZMIENNE GLOBALNE --------------#

m = 1
E = -0.6
v0 = 0

#-------------------- FUNKCJE GLOBALNE --------------#

#Pochodna pierwszego rzedu
def diff(x,dx,fun):
	return (fun(x+dx)-fun(x-dx))/(2*dx)

#Funkcja potencjalu V(x)
def potencjal(x):
	return -numpy.exp(-x**2) -1.2*numpy.exp(-(x-2)**2)

#Rysowanie potencjalu V(x)

def rysuj_potencjal(a,b,fun,fi,name):
	if fi:
		x = numpy.linspace(a,b,1000)
		fig, ax = plt.subplots()
		ax.plot(x,fun(x))
		ax.grid()
		ax.set(xlabel = 'x', ylabel = 'V(x)')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

#-------------------- ZADANIE 1 ---------------------#

#Funkcja F
def funkcja_f1(x,E_0 = E):
	return potencjal(x)-E_0

rysuj_potencjal(-2,4,potencjal,draw_chart,'pot')

#____________________ ZADANIE 1.1 ___________________#


#Funckja liczaca miejsce zerowe funkcji fun metoda bisekcji
def bisekcja(err,a,b,fun):
	if not ((fun(a)>0 and fun(b)<0) or (fun(a)<0 and fun(b)>0)):
		print('Podaj inne wartości graniczne:', a, b)
	else:
		TAB_x = [a,b]
		TAB_s = ['a','b']
		n = 0
		while(abs(a - b) > err):
			mid = (a + b)/2
			n += 1
			TAB_x.append(mid)
			TAB_s.append(str(n))
			if fun(a) * fun(mid) < 0:
				b = mid
			elif fun(b) * fun(mid) < 0:
				a = mid

		#Zwrot wynikow
		zad = '#____________________ ZADANIE 1.1 ___________________#'
		head = 'a = ' + "{: .6f}".format(TAB_x[0]) + ' b = ' + "{: .6f}".format(TAB_x[1]) + ' ERR = ' + "{: .6f}".format(err)
		print(zad, '\n', '\n', head)
		file.write(zad + '\n' + '\n' + head + '\n')
		for i in list(zip(TAB_x, TAB_s)):
			out = '{: .6f}'.format(i[0]) + '\t' + "{: .6f}".format(fun(i[0])) + '\t' + i[1]
			print(out)
			file.write(out + '\n')
		return [TAB_x,TAB_s,n]

#Funkcje do rysowania wykresów
def rysuj_1_1(wynik,a,b,fi,legend,name):
	if fi:
		s = wynik[1]
		x = numpy.array(wynik[0])
		xx = numpy.linspace(a,b,1000)
		fig, ax = plt.subplots()
		ax.plot(xx, funkcja_f1(xx), color = 'green')
		ax.plot(x,funkcja_f1(x), 'rd',)
		for i in zip(s,x):
			ax.text(i[1],0.03 + funkcja_f1(i[1]), i[0])
		ax.grid()
		ax.set(xlabel = 'x [m]', ylabel = 'F(x) [J]')
		punkty = mlines.Line2D([], [], label = 'Punkty z metody bisekcji', color = 'red', linestyle = 'None', marker = 'd')
		linia = mlines.Line2D([], [], label = 'Funkcja F(x)', color = 'green', linestyle = 'solid', marker = 'None')
		ax.legend(handles = [punkty, linia], loc = legend)
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

def rysuj_1_2(wynik,fi,name):
	if fi:
		x = numpy.array(wynik[0])
		n = numpy.arange(1,len(x)+1,1)
		fig, ax = plt.subplots()
		ax.semilogy(n,abs(funkcja_f1(x)), 'rd',)
		ax.grid()
		ax.set(xlabel = 'Liczba iteracji i', ylabel = '|F(x)|')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

#Punkty poczatkowe
a_1 = 1
b_1 = 3
a_2 = -2
b_2 = 0

#Dopuszczalny blad err
err = 1e-6

#Wynik 1
WYNIK_1_1 = bisekcja(err,a_1,b_1,funkcja_f1)
#Wynik 2
WYNIK_1_2 = bisekcja(err,a_2,b_2,funkcja_f1)

#Rysowanie wykresow
rysuj_1_1(WYNIK_1_1, 0, 4, draw_chart,'upper left','1_1_1_1')
rysuj_1_2(WYNIK_1_1,draw_chart, '1_1_1_2')

rysuj_1_1(WYNIK_1_2, -3, 0, draw_chart,'lower left', '1_1_2_1')
rysuj_1_2(WYNIK_1_2,draw_chart, '1_1_2_2')



#____________________ ZADANIE 1.2 ___________________#


#Funckja liczaca miejsce zerowe funkcji fun metoda Newtona-Raphsona
def new_rap(err,x_n,fun,dx):
	x_n1 = x_n - fun(x_n)/diff(x_n,dx,fun)
	TAB = [x_n, x_n1]
	
	while(abs(x_n - x_n1) > err):
		x_n = x_n1
		x_n1 = x_n - fun(x_n)/diff(x_n,dx,fun)
		TAB.append(x_n1)

	#Zwrot wynikow
	zad = '#____________________ ZADANIE 1.2 ___________________#'
	head = 'x0 = ' + "{: .6f}".format(TAB[0]) + ' ERR = ' + "{: .6f}".format(err)
	print(zad, '\n', '\n', head)
	file.write(zad + '\n' + '\n' + head + '\n')
	for i in TAB:
		out = '{: .6f}'.format(i) + '\t' + "{: .6f}".format(fun(i))
		print(out)
		file.write(out + '\n')
	return TAB

#Funckje do rysowania wykresów
def rysuj_1_11(wynik,a,b,fi,legend,name):
	if fi:
		s = numpy.arange(0,len(wynik),1)
		x = numpy.array(wynik)
		xx = numpy.linspace(a,b,1000)
		fig, ax = plt.subplots()
		ax.plot(xx, funkcja_f1(xx), color = 'green')
		ax.plot(x,funkcja_f1(x), 'rd',)
		for i in zip(s,x):
			ax.text(i[1],0.03 + funkcja_f1(i[1]), i[0])
		ax.grid()
		ax.set(xlabel = 'x [m]', ylabel = 'F(x) [J]')
		punkty = mlines.Line2D([], [], label = 'Punkty z metody Newtona-Raphsona', color = 'red', linestyle = 'None', marker = 'd')
		linia = mlines.Line2D([], [], label = 'Funkcja F(x)', color = 'green', linestyle = 'solid', marker = 'None')
		ax.legend(handles = [punkty, linia], loc = legend)
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

def rysuj_1_22(wynik,fi,name):
	if fi:
		x = numpy.array(wynik)
		n = numpy.arange(1,len(x)+1,1)
		fig, ax = plt.subplots()
		ax.semilogy(n,abs(funkcja_f1(x)), 'rd',)
		ax.grid()
		ax.set(xlabel = 'Liczba iteracji i', ylabel = '|F(x)|')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()


#Punkty poczatkowe
a = 3
b = 1
dx = 0.01
#Dopuszczalny blad err
err = 1e-6

WYNIK_1_1 = new_rap(err,a,funkcja_f1,dx)
WYNIK_1_2 = new_rap(err,b,funkcja_f1,dx)

#Rysowanie wykresow
rysuj_1_11(WYNIK_1_1, 0, 4, draw_chart,'upper left','1_2_1_1')
rysuj_1_22(WYNIK_1_1,draw_chart, '1_2_1_2')

rysuj_1_11(WYNIK_1_2, -3, 0, draw_chart,'lower left', '1_2_2_1')
rysuj_1_22(WYNIK_1_2,draw_chart, '1_2_2_2')



#-------------------- ZADANIE 2 ---------------------#

#Funkcje do jawnego schematu Eulera

def predkosc(v0, x0, dt, dx, alpha = 0):
	return v0 - 1/m * diff(x0,dx,potencjal)*dt - v0*alpha*dt

def polozenie(x0, v0, dt):
	return x0 + v0 * dt

def energia_kinetyczna(v):
	def energia(v):
		return m*v**2/2
	return list(map(energia,v))

def energia_potencjalna(x):
	return list(map(potencjal,x))

def energia_cal(Ek, V):
	return list(map(lambda x,y : x+y, Ek, V))

def funkcja_f2(t0, tk, v0, x0, dt, dx, alpha = 0):
	t = numpy.arange(t0, tk, dt)
	x = [x0]
	v = [v0]
	for i in range(len(t)-1):
		vn = predkosc(v0,x0,dt,dx,alpha)
		v.append(vn)
		v0 = vn
		xn = polozenie(x0,v0,dt)
		x.append(xn)
		x0 = xn

	return (x,v,t)

def rysuj_2(x,y,fi,name,x_lab, y_lab):
	if fi:
		fig, ax = plt.subplots()
		ax.plot(x,y, color = 'blue')

		ax.grid()
		ax.set(xlabel = x_lab, ylabel = y_lab)
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()



#Punkty poczatkowe 2_1
t0 = 0
tk = 30
x0 = WYNIK_1_1[-1]
dt = 0.01
dx = 0.01
#Wyniki 2_1
WYNIK_2_1 = funkcja_f2(t0, tk, v0, x0, dt, dx)
rysuj_2(WYNIK_2_1[2], WYNIK_2_1[0], draw_chart, '2_1x', 't [s]', 'x [m]')
rysuj_2(WYNIK_2_1[2], WYNIK_2_1[1], draw_chart, '2_1v', 't [s]', 'v [m/s]')
E_k = energia_kinetyczna(WYNIK_2_1[1])
rysuj_2(WYNIK_2_1[2], E_k, draw_chart, '2_1ek', 't [s]', '$E_k(t)$ [J]')
V = energia_potencjalna(WYNIK_2_1[0])
rysuj_2(WYNIK_2_1[0], V, draw_chart, '2_1vv', 'x [m]', '$V(x(t))$ [J]')
Ec = energia_cal(E_k,V)
rysuj_2(WYNIK_2_1[2], Ec, draw_chart, '2_1ekv', 't [s]', '$E_k(t) + V(t)$ [J]')

#Punkty poczatkowe 2_2_1
t0 = 0
tk = 100
x0 = WYNIK_1_1[-1]
dt = 0.01
dx = 0.01
#Wyniki 2_2_1
WYNIK_2_2_1 = funkcja_f2(t0, tk, v0, x0, dt, dx)
rysuj_2(WYNIK_2_2_1[0], WYNIK_2_2_1[1], draw_chart, '2_2_1', 'x [m]', 'v [m/s]')

#Punkty poczatkowe 2_2_2
t0 = 0
tk = 1000
x0 = WYNIK_1_1[-1]
dt = 0.001
dx = 0.01
#Wyniki 2_2_2
WYNIK_2_2_2 = funkcja_f2(t0, tk, v0, x0, dt, dx)
rysuj_2(WYNIK_2_2_2[0], WYNIK_2_2_2[1], draw_chart, '2_2_2', 'x [m]', 'v [m/s]')



#-------------------- ZADANIE 3 ---------------------#

def rysuj_3(x,y,fi,name,x_lab, y_lab):
	if fi:
		fig, ax = plt.subplots()
		ax.semilogy(x,y, color = 'blue')

		ax.grid()
		ax.set(xlabel = x_lab, ylabel = y_lab)
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

#Punkty poczatkowe 3_3
t0 = 0
tk = 30
x0 = WYNIK_1_1[-1]
dt = 0.01
dx = 0.01

#Punkty poczatkowe 3_1
alpha = 0.5
#Wyniki 3_1
WYNIK_3_1 = funkcja_f2(t0, tk, v0, x0, dt, dx, alpha)
rysuj_2(WYNIK_3_1[0], WYNIK_3_1[1], draw_chart, '3_1', 'x [m]', 'v [m/s]')

#Punkty poczatkowe 3_2
alpha = 5
#Wyniki 3_2
WYNIK_3_2 = funkcja_f2(t0, tk, v0, x0, dt, dx, alpha)
rysuj_2(WYNIK_3_2[0], WYNIK_3_2[1], draw_chart, '3_2', 'x [m]', 'v [m/s]')

#Punkty poczatkowe 3_3
alpha = 201
#Wyniki 3_3
WYNIK_3_3 = funkcja_f2(t0, tk, v0, x0, dt, dx, alpha)
rysuj_3(WYNIK_3_3[0], WYNIK_3_3[1], draw_chart, '3_3', 'x [m]', 'v [m/s]')



