#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



#Rysowanie wykresów
draw_chart = (len(sys.argv)>1)

#-------------------- ZMIENNE GLOBALNE --------------#

au = 149597870700 	#m
m = 1.989e30		#kg
G = 6.6741e-11		#m^3/kg*s^2
v_x0 = 54600		#m/s
v_y0 = 0			#m/s
x_0 = 0 * au		#au -> m
y_0 = 0.586 * au	#au -> m
T = 100				#y

#-------------------- FUNKCJE GLOBALNE --------------#

#Przyspieszenie a_x
def a_x(x,y, G = G):
	r = numpy.sqrt(x**2 + y**2)
	return -G*m/(r**3)*x

#Przyspieszenie a_y
def a_y(x,y, G = G):
	r = numpy.sqrt(x**2 + y**2)
	return -G*m/(r**3)*y

#Funkcja licząca ilosc krokow czasowych: T -> lata
def N_krok(T, n, dt):
	return int(T*365*24*3600*n/dt) + 1

#Funkcje do rysowania wykresów
def rysuj_yx(wynik, name, fi = draw_chart):
	if fi:
		X = numpy.array(wynik[0])/au
		Y = numpy.array(wynik[1])/au
		fig, ax = plt.subplots()
		ax.plot(X,Y)
		ax.grid()
		ax.set(xlabel = 'x [au]', ylabel = 'y [au]')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

def rysuj_yt(wynik, name, fi = draw_chart):
	if fi:
		Y = numpy.array(wynik[1])/au
		t = numpy.array(wynik[2])/(365*24*3600)
		fig, ax = plt.subplots()
		ax.plot(t,Y)
		ax.grid()
		ax.set(xlabel = 't [year]', ylabel = 'y [au]')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()

def rysuj_ydt(wynik, name, fi = draw_chart):
	if fi:
		X = numpy.array(wynik[0])/au
		Y = numpy.array(wynik[1])/au
		R = numpy.sqrt(X**2 + Y**2)
		dt = numpy.array(wynik[3])/(3600)
		fig, ax = plt.subplots()
		ax.semilogy(R[1:], dt)
		ax.grid()
		ax.set(ylabel = '$\Delta$ t [h]', xlabel = 'R [au]')
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)
		fig.savefig(name+'.png', dpi = 100)
		plt.show()


#-------------------- ZADANIE 1 ---------------------#

#Funkcja liczaca oribite komety metoda jawnego schematu Eulera
def jawny_euler(dt, x_0 = x_0, y_0 = y_0, v_x0 = v_x0, v_y0 = v_y0, n = 3):

	def polozenie(r, v, dt):
		return r + v*dt

	def predkosc(v, a, dt):
		return v + a*dt
	t = [0]
	X = [x_0]
	Y = [y_0]
	v_x = v_x0
	v_y = v_y0

	for i in range(N_krok(T,n,dt)):
		x = X[i] 
		y = Y[i]
		v_yn = predkosc(v_y, a_y(x, y), dt)
		v_xn = predkosc(v_x, a_x(x, y), dt)
		y_n = polozenie(y, v_y, dt)
		x_n = polozenie(x, v_x, dt)

		v_y = v_yn
		v_x = v_xn
		Y.append(y_n)
		X.append(x_n)
		t.append(t[i]+dt)

	return (X,Y,t)



def zad_1():
	dt = 60*60	#60 min
	t0 = time.time()
	R = jawny_euler(dt)
	print('#-------------------- ZADANIE 1 ---------------------# \n')
	print('TIME [s] :	', "{: .6f}".format(time.time() - t0),'\n')
	rysuj_yx(R, "1_yx")
	rysuj_yt(R, "1_yt")

#zad_1()


#-------------------- ZADANIE 2 ---------------------#

def RK4(dt, x_0 = x_0, y_0 = y_0, v_x0 = v_x0, v_y0 = v_y0, n = 3):
	
	def fun(v_x , v_y, x , y , dt = dt):

		def fin(u, k1, k2, k3, k4, dt = dt):
			return u + dt/6 * (k1 + 2*k2 + 2*k3 + k4) 

		k_11 = v_x
		k_12 = v_y
		k_13 = a_x(x,y)
		k_14 = a_y(x,y)


		k_21 = v_x + dt/2 * k_13
		k_22 = v_y + dt/2 * k_14
		k_23 = a_x(x + dt/2 * k_11, y + dt/2 * k_12)
		k_24 = a_y(x + dt/2 * k_11, y + dt/2 * k_12)

		k_31 = v_x + dt/2 * k_23
		k_32 = v_y + dt/2 * k_24
		k_33 = a_x(x + dt/2 * k_21, y + dt/2 * k_22)
		k_34 = a_y(x + dt/2 * k_21, y + dt/2 * k_22)

		k_41 = v_x + dt * k_33
		k_42 = v_y + dt * k_34
		k_43 = a_x(x + dt * k_31, y + dt * k_32)
		k_44 = a_y(x + dt * k_31, y + dt * k_32)	

		v_xn = fin(v_x, k_13, k_23, k_33, k_43)
		v_yn = fin(v_y, k_14, k_24, k_34, k_44)
		x_n = fin(x, k_11, k_21, k_31, k_41)
		y_n = fin(y, k_12, k_22, k_32, k_42)

		return (x_n, y_n, v_xn, v_yn)

	t = [0]
	X = [x_0]
	Y = [y_0]
	v_x = v_x0
	v_y = v_y0

	for i in range(N_krok(T,n,dt)):
		x = X[i]
		y = Y[i]

		wynik = fun(v_x, v_y, x, y)

		v_y = wynik[3]
		v_x = wynik[2]
		y_n = wynik[1]
		x_n = wynik[0]
		Y.append(y_n)
		X.append(x_n)
		t.append(t[i]+dt)

	return (X,Y,t)

def zad_2():
	dt = 60*60	#60 min
	t0 = time.time()
	R = RK4(dt)
	print('#-------------------- ZADANIE 2---------------------# \n')
	print('TIME [s] :	', "{: .6f}".format(time.time() - t0),'\n')
	rysuj_yx(R, "2_yx")
	rysuj_yt(R, "2_yt")

#zad_2()



#-------------------- ZADANIE 3 ---------------------#

#Funkcja liczaca oribite komety metoda jawnego schematu Eulera z automatycznym doborem kroku czasowego
def jawny_euler_at(dt_0, tol, x_0 = x_0, y_0 = y_0, v_x0 = v_x0, v_y0 = v_y0, n = 3):

	def h_eps(a, b, N =1):
		return (b - a)/(2**N -1)

	def polozenie(r, v, dt):
		return r + v*dt

	def predkosc(v, a, dt):
		return v + a*dt

	t = [0]
	tt = 0
	D = []
	X = [x_0]
	Y = [y_0]
	v_x = v_x0
	v_y = v_y0
	dt = dt_0
	epsilon = 0
	i = 0
	c = 0.9
	N = 1
	while 3*T*365*24*3600 >= tt:
		x = X[i]
		y = Y[i]
		while abs(epsilon) > tol or epsilon == 0:

			y_n = polozenie(y, v_y, dt)
			x_n = polozenie(x, v_x, dt)

			y_nn = polozenie(y, v_y, dt/2)
			x_nn = polozenie(x, v_x, dt/2)

			v_ynn = predkosc(v_y, a_y(x_nn, y_nn), dt/2)
			v_xnn = predkosc(v_x, a_x(x_nn, y_nn), dt/2)

			y_nnn = polozenie(y_nn, v_ynn, dt/2)
			x_nnn = polozenie(x_nn, v_xnn, dt/2)

			epsilon_x = h_eps(x_n, x_nnn)
			epsilon_y = h_eps(y_n, y_nnn)

			epsilon = epsilon_x if abs(epsilon_x) > abs(epsilon_y) else epsilon_y
			dt_n = c*dt * abs((tol/epsilon))**(1/2)
			dt = dt_n


		v_yn = predkosc(v_y, a_y(x, y), dt)
		v_xn = predkosc(v_x, a_x(x, y), dt)
		y_n = polozenie(y, v_y, dt)
		x_n = polozenie(x, v_x, dt)

		v_y = v_yn
		v_x = v_xn

		Y.append(y_n)
		X.append(x_n)
		t.append(t[i]+dt)
		D.append(dt)
		tt += dt
		i += 1
		#Restart pętli
		epsilon = 0
	return (X,Y,t,D)



def zad_3(tol):
	dt = 60 * 60	#60 min
	t0 = time.time()
	R = jawny_euler_at(dt, tol)
	print('#-------------------- ZADANIE 3 ---------------------# \n')
	print('TIME [s] :	', "{: .6f}".format(time.time() - t0),'\n')
	rysuj_yx(R, "3_yx_" + str(tol))
	rysuj_yt(R, "3_yt_" + str(tol))
	rysuj_ydt(R, "3_ydt_" + str(tol))

#zad_3(1000)
#zad_3(100)



#-------------------- ZADANIE 4 ---------------------#

def RK4_at(dt_0, tol, x_0 = x_0, y_0 = y_0, v_x0 = v_x0, v_y0 = v_y0, n = 3):
	
	dt = dt_0

	def h_eps(a, b, N = 4):
		return (b - a)/((2**N) -1)

	def fun(v_x , v_y, x , y , dt = dt):

		def fin(u, k1, k2, k3, k4, dt = dt):
			return u + dt/6 * (k1 + 2*k2 + 2*k3 + k4) 

		k_11 = v_x
		k_12 = v_y
		k_13 = a_x(x,y)
		k_14 = a_y(x,y)


		k_21 = v_x + dt/2 * k_13
		k_22 = v_y + dt/2 * k_14
		k_23 = a_x(x + dt/2 * k_11, y + dt/2 * k_12)
		k_24 = a_y(x + dt/2 * k_11, y + dt/2 * k_12)

		k_31 = v_x + dt/2 * k_23
		k_32 = v_y + dt/2 * k_24
		k_33 = a_x(x + dt/2 * k_21, y + dt/2 * k_22)
		k_34 = a_y(x + dt/2 * k_21, y + dt/2 * k_22)

		k_41 = v_x + dt * k_33
		k_42 = v_y + dt * k_34
		k_43 = a_x(x + dt * k_31, y + dt * k_32)
		k_44 = a_y(x + dt * k_31, y + dt * k_32)	

		v_xn = fin(v_x, k_13, k_23, k_33, k_43)
		v_yn = fin(v_y, k_14, k_24, k_34, k_44)
		x_n = fin(x, k_11, k_21, k_31, k_41)
		y_n = fin(y, k_12, k_22, k_32, k_42)

		return (x_n, y_n, v_xn, v_yn)

	t = [0]
	tt = 0
	D = []
	X = [x_0]
	Y = [y_0]
	v_x = v_x0
	v_y = v_y0

	epsilon = 0
	i = 0
	c = 0.9
	N = 4

	while n*T*365*24*3600 >= tt:
		x = X[i]
		y = Y[i]
		while abs(epsilon) > tol or epsilon == 0:
			wynik = fun(v_x, v_y, x, y, dt)
			y_n = wynik[1]
			x_n = wynik[0]

			wynik_nn = fun(v_x, v_y, x, y, dt/2)
			x_nn = wynik_nn[0]
			y_nn = wynik_nn[1]
			v_xnn = wynik_nn[2]
			v_ynn = wynik_nn[3]

			wynik_nnn = fun(v_xnn, v_ynn, x_nn, y_nn, dt/2)
			y_nnn = wynik_nnn[1]
			x_nnn = wynik_nnn[0]		

			epsilon_x = h_eps(x_n, x_nnn)
			epsilon_y = h_eps(y_n, y_nnn)

			epsilon = epsilon_x if abs(epsilon_x) > abs(epsilon_y) else epsilon_y
			dt_n = c*dt * abs((tol/epsilon))**(1/(N + 1))
			dt = dt_n
		
		wynik = fun(v_x, v_y, x, y, dt)

		v_y = wynik[3]
		v_x = wynik[2]
		y_n = wynik[1]
		x_n = wynik[0]

		Y.append(y_n)
		X.append(x_n)
		t.append(t[i]+dt)
		D.append(dt)
		tt += dt
		i += 1
		#Restart pętli
		epsilon = 0
	return (X,Y,t,D)



def zad_4(tol):
	dt = 60 * 60	#60 min
	t0 = time.time()
	R = RK4_at(dt, tol)
	print('#-------------------- ZADANIE 4 ---------------------# \n')
	print('TIME [s] :	', "{: .6f}".format(time.time() - t0),'\n')
	rysuj_yx(R, "4_yx_" + str(tol))
	rysuj_yt(R, "4_yt_" + str(tol))
	rysuj_ydt(R, "4_ydt_" + str(tol))


#zad_4(1000)
#zad_4(100)
#zad_4(1)












