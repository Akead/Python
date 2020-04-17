#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



#Rysowanie wykresów
draw_chart = (len(sys.argv)>1)

#-------------------- ZMIENNE GLOBALNE --------------#

m = 1		#kg
g = 9.81	#m/s^2
R = 1		#m
dt = 0.01 	#s
N = 1000
th_dt_0 = 0 #rad/s

#-------------------- FUNKCJE GLOBALNE --------------#

#Rysowanie wykresów
def rysuj(lista, arg_1, arg_2, name_x, name_y, name, fi = draw_chart):
	fig, ax = plt.subplots()
	legend = []
	for i in lista:
		A = i.wyp()
		S = i.chart_par()
		ax.plot(A[arg_1], A[arg_2], color = S[1], linestyle = S[2])
		opis = mlines.Line2D([], [], label = S[0], linestyle = S[2], color = S[1])
		legend.append(opis)
	ax.legend(handles = legend, loc = 'upper right', fontsize = 'x-large')
	ax.set(ylabel = name_y, xlabel = name_x)

	ax.grid()
	ax.xaxis.label.set_size(18)
	ax.yaxis.label.set_size(18)
	fig = plt.gcf()
	fig.set_size_inches(18.5, 10.5, forward = True)
	fig.savefig(name+'.png', dpi = 100)

	if fi:
		plt.show()



#Zamiana stopni na radiany
def deg_to_rad(x):
	return x/360*2*numpy.pi


#Przyspieszenie katowe
def theta_2dt(theta, R = R, g = g):
	return -g/R * numpy.sin(theta)

#Metoda RK4
def RK4(th_0, th_dt_0 = th_dt_0, dt = dt, N = N):

	def fun(th, th_dt, dt = dt):

		def fin(u, k1, k2, k3, k4, dt = dt):
			return u + dt/6 * (k1 + 2*k2 + 2*k3 + k4) 

		k_11 = th_dt
		k_13 = theta_2dt(th)

		k_21 = th_dt + dt/2 * k_13
		k_23 = theta_2dt(th + dt/2 * k_11)

		k_31 = th_dt + dt/2 * k_23
		k_33 = theta_2dt(th + dt/2 * k_21)

		k_41 = th_dt + dt * k_33
		k_43 = theta_2dt(th + dt * k_31)


		th_dt_n = fin(th_dt, k_13, k_23, k_33, k_43)

		th_n = fin(th, k_11, k_21, k_31, k_41)


		return (th_n, th_dt_n)

	t = [0]
	Th = [th_0]
	Th_dt = [th_dt_0]

	for i in range(N):
		th = Th[i]
		th_dt = Th_dt[i]

		wynik = fun(th, th_dt)
		th_n = wynik[0]
		th_dt_n = wynik[1]
		Th.append(th_n)
		Th_dt.append(th_dt_n)
		t.append(t[i]+dt)

	return (Th, Th_dt ,t)

#-------------------- KLASY --------------#

class wahadlo:
	#Wahadło i jego parametry

	Th = []
	Th_dt = []
	t = []
	T = []
	U = []
	E = []
	Okres = 0
	

	def __init__(self, nazwa, kolor, theta_0, linetype = 'solid' ,theta_dt_0 = th_dt_0, N = N, dt = dt, R = R, m = m):
		self.nazwa = nazwa
		self.kolor = kolor
		self.th_0 = theta_0
		self.th_dt_0 = theta_dt_0
		self.N = N
		self.dt = dt
		self.R = R
		self.m = m
		self.linetype = linetype

	def cal(self):
		self.Th, self.Th_dt, self.t = RK4(self.th_0, self.th_dt_0 , self.dt, self.N)
	
	def cycle(self):
		s = self.Th[0]
		for i in self.Th:
			self.Okres += self.dt
			if self.Okres > 1 and s-i < 1e-3:
				break


	def analitycznie(self):
		omega = numpy.sqrt(g/self.R)
		t = 0
		for i in range(self.N):
			self.t.append(dt*i)

		def h_fun(time):
			return self.th_0 * numpy.sin(omega * time + numpy.pi/2)
		self.Th = list(map(h_fun, self.t))

	def oblicz(self):
		def kinetyczna(phi_dt):
			return self.m/2 * self.R**2 * phi_dt**2

		def potencjalna(phi):
	 		return - self.m * g * self.R * numpy.cos(phi)

		self.T = list(map(kinetyczna, self.Th_dt))
		self.U = list(map(potencjalna, self.Th))
		self.E = list(map(lambda x,y : x+y, self.T, self.U))

	def wyp(self):
		return (numpy.array(self.Th), numpy.array(self.t), numpy.array(self.Th_dt), numpy.array(self.T), numpy.array(self.U), numpy.array(self.E), self.Okres, self.th_0)

	def chart_par(self):
		return (self.nazwa, self.kolor, self.linetype)



#-------------------- ZADANIE 1 ---------------------#
def zadanie_1():
	A = wahadlo("Numerycznie: $\phi(t=0)$ = 4°", "red", deg_to_rad(4), (0,(5,10)))
	A.cal()
	B = wahadlo("Analitycznie: $\phi(t=0)$ = 4°", "green", deg_to_rad(4))
	B.analitycznie()

	l = [A,B]

	rysuj(l,1,0,"t [s]", '$\phi$ [rad]', "zad_1")

#zadanie_1()

#-------------------- ZADANIE 2 ---------------------#

def zadanie_2():

	A4 = wahadlo("Numerycznie: $\phi(t=0)$ = 4°", "red", deg_to_rad(4))
	A45 = wahadlo("Numerycznie: $\phi(t=0)$ = 45°", "green", deg_to_rad(45))
	A90 = wahadlo("Numerycznie: $\phi(t=0)$ = 90°", "blue", deg_to_rad(90))
	A135 = wahadlo("Numerycznie: $\phi(t=0)$ = 135°", "orange", deg_to_rad(135))
	A175 = wahadlo("Numerycznie: $\phi(t=0)$ = 175°", "purple", deg_to_rad(175))

	L = [A4, A45, A90, A135, A175]

	for i in L:
		i.cal()
		i.oblicz()



	#phi(t)
	rysuj(L,1,0,"t [s]", '$\phi$ [rad]', "zad_2_phi_t")

	#T(t)
	rysuj(L,1,3,"t [s]", 'T [J]', "zad_2_T")

	#U(t)
	rysuj(L,1,4,"t [s]", 'U [J]', "zad_2_U")

	#E(t)
	rysuj(L,1,5,"t [s]", 'E [J]', "zad_2_E")

	#Z(t)
	rysuj(L,0,2,"$\phi$ [rad]", '$\dot{\phi} \, [\\frac{rad}{s}]$', "zad_2_Z")

#zadanie_2()


#-------------------- ZADANIE 3 ---------------------#

def zadanie_3():

	data = numpy.linspace(1e-3,numpy.pi-1e-3,360)
	L = [wahadlo('','', i) for i in data]
	
	phi = []
	T = []

	for i in L:
		i.cal()
		i.cycle()
		Q = i.wyp()
		T.append(Q[-2])
		phi.append(Q[-1])
	phi = numpy.array(phi)
	T = numpy.array(T)

	fig, ax = plt.subplots()
	legend = []
	ax.plot(phi, T, color = 'red')
	opis = mlines.Line2D([], [], label = 'Uzyskane dane',color = 'red')
	legend.append(opis)
	ax.legend(handles = legend, loc = 'upper left', fontsize = 'x-large')
	ax.set(ylabel = 'T [s]', xlabel = '$\phi$ [rad]')

	ax.grid()
	ax.xaxis.label.set_size(18)
	ax.yaxis.label.set_size(18)
	fig = plt.gcf()
	fig.set_size_inches(18.5, 10.5, forward = True)
	fig.savefig('zadanie_3.png', dpi = 100)

	if draw_chart:
		plt.show()

#zadanie_3()











