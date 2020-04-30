#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



#Chart drawing
draw_chart = (len(sys.argv)<=1)


#-------------------- GLOBAL VARIABLES --------------#

N = 5000
q = B = m = 1



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

class charged_particle:

	def __init__(self, N, start_cond, name = '', color = 'blue', marker = None, line_type = 'solid', linewidth = 1, q = q, B = B, m = m):
		self.__N = N
		self.__name = name
		self.__color = color
		self.__marker = marker
		self.__line_type = line_type
		self.__linewidth = linewidth
		self.__q = q
		self.__B = B
		self.__m = m
		self.__omega_c = q*B/m
		self.__T = 2*numpy.pi/self.__omega_c
		self.__dt = 5*self.__T/self.__N

		self.__resoult = [start_cond]
		self.__t = [0]
		self.__Energy = numpy.array([])

	def __function(self, s_vector):
		s_0 = s_vector[0]
		s_1 = s_vector[1]
		s_2 = s_vector[2]
		s_3 = s_vector[3]
		s_4 = s_vector[4]
		s_5 = s_vector[5]

		d_0 = s_3/self.__m
		d_1 = s_4/(self.__m*s_0**2) - self.__q*self.__B/(2*self.__m)
		d_2 = s_5/self.__m
		d_3 = (s_4**2)/(self.__m*s_0**3) - self.__q**2 * self.__B**2 *s_0/(4*self.__m)
		d_4 = 0
		d_5 = 0

		return numpy.array([d_0, d_1, d_2, d_3, d_4, d_5])

	def __cylinder_to_kartesian(self):
		tmp = numpy.rot90(self.__resoult)
		r = tmp[5]
		phi = tmp[4]
		z = tmp[3]


		x = r * numpy.cos(phi)
		y = r * numpy.sin(phi)
		z = z

		return (x,y,z)

	def __total_energy(self):
		tmp = numpy.rot90(self.__resoult)
		r = tmp[5]
		phi = tmp[4]
		z = tmp[3]
		p_r = tmp[2]
		p_f = tmp[1]
		p_z = tmp[0]

		A = 1/(2*self.__m) * (p_r**2 + p_f**2/r**2 + p_z**2)
		B = self.__q * self.__B/(2*self.__m) * p_f
		C = self.__q**2 * self.__B**2 * r**2/(8*self.__m)
		self.__Energy = A - B + C

	def __plot(self, x, y, x_label, y_label, legend_label, name, color = None, marker = None, linetype = None, draw = draw_chart, save = True):
		if color == None:
			color = self.__color
		if marker == None:
			marker = self.__marker
		if linetype == None:
			linetype = self.__line_type
	
		fig, ax = plt.subplots()
		ax.plot(x,y, color = color, marker = marker, linestyle = linetype)

		#Legend
		legend = []
		description = mlines.Line2D([], [], label = legend_label, color = color, marker = marker, linestyle = linetype)
		legend.append(description)
		ax.legend(handles = legend, loc = 'lower right', fontsize = 'x-large')

		#Label XY
		ax.set(ylabel = y_label, xlabel = x_label)

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if save:
			fig.savefig(name+'.png', dpi = 100)

		if draw:
			plt.show()

	def calculate(self):
		
		for i in range(self.__N):
			start_conditions = self.__resoult[i]
			new_conditions = RK4(start_conditions, self.__function, self.__dt)
			self.__resoult.append(new_conditions)
			self.__t.append(self.__t[i] + self.__dt)

		self.__total_energy()

	def return_resoults(self):
		tmp = numpy.rot90(self.__resoult)
		r = tmp[5]
		phi = tmp[4]
		z = tmp[3]
		p_r = tmp[2]
		p_f = tmp[1]
		p_z = tmp[0]
		x, y, z = self.__cylinder_to_kartesian()

		return [(r, phi, z),(p_r, p_f, p_z),(x, y, z), self.__Energy, self.__t]

	def return_properities(self):
		return (self.__name, self.__color, self.__marker, self.__line_type, self.__linewidth)

	def print_resoults(self):

		print(80*'-')
		print('r     ', '	', '	phi     ','	', 'z   ', '	','\t p_z   ', '	','p_phi   ','	','p_z   ','	', 'E')	
		print(80*'-')
		for i, k in zip(self.__resoult, self.__Energy):
			print('{:.6f}'.format(i[0]), '	','{:.6f}'.format(i[1]), '	', '{:.6f}'.format(i[2]), '	', '{:.6f}'.format(i[3]),  '	', '{:.6f}'.format(i[4]), '	', '{:.6f}'.format(i[5]),'	', '{:.6f}'.format(k), '\n')
		print(80*'-')

	def plot(self):
	
		tmp = numpy.rot90(self.__resoult)
		r = tmp[5]
		phi = tmp[4]
		z = tmp[3]
		p_r = tmp[2]
		p_f = tmp[1]
		p_z = tmp[0]

		x, y, z = self.__cylinder_to_kartesian()

		self.__plot(x, y, 'x(t)', 'y(t)', '$\phi (t)$', self.__name + 'xy_t')
		self.__plot(self.__t, self.__Energy, 't', 'E', 'E(t)', self.__name + 'E_t')
		self.__plot(self.__t, r, 't', 'r', 'r(t)', self.__name + 'r_t')
		self.__plot(self.__t, phi, 't', '$\phi$', '$\phi (t)$', self.__name + 'phi_t')
		self.__plot(self.__t, p_r, 't', '$p_r$', '$\omega (t)$', self.__name + 'p_r__t')
		self.__plot(self.__t, p_f, 't', '$p_{\phi}$', '$\omega (t)$', self.__name + 'p_phi__t')


class draw_plots:

	def __init__(self, L):
		self.__L = L
		self.__names = [i.return_properities()[0] for i in self.__L]
		self.__colors = [i.return_properities()[1] for i in self.__L]
		self.__markers = [i.return_properities()[2] for i in self.__L]
		self.__lines = [i.return_properities()[3] for i in self.__L]
		self.__linewidth = [i.return_properities()[4] for i in self.__L]

		self.__cylinder = [i.return_resoults()[0] for i in self.__L]
		self.__momentums = [i.return_resoults()[1] for i in self.__L]
		self.__kartesian = [i.return_resoults()[2] for i in self.__L]
		self.__Energy = [i.return_resoults()[3] for i in self.__L]
		self.__t = [i.return_resoults()[4] for i in self.__L]

	def __plot(self, X, Y, x_label, y_label, label, save_name, draw = draw_chart, save = True):

		#Figure
		fig, ax = plt.subplots()

		#Legend
		legend = []

		for x, y, name, color, marker, linetype, linewidth in zip(X, Y, self.__names, self.__colors, self.__markers, self.__lines, self.__linewidth):
			ax.plot(x,y, color = color, marker = marker, linestyle = linetype, linewidth = linewidth, markersize = 5)
			description = mlines.Line2D([], [], label = name + ' : ' + label, color = color, marker = marker, linestyle = linetype)
			legend.append(description)

		ax.legend(handles = legend, loc = 'lower right', fontsize = 'x-large')

		#Label XY
		ax.set(ylabel = y_label, xlabel = x_label)

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		fig = plt.gcf()
		fig.set_size_inches(18.5, 10.5, forward = True)

		if save:
			fig.savefig(save_name+'.png', dpi = 100)

		if draw:
			plt.show()

	def plot(self):
		x = numpy.rot90(self.__kartesian)[2]
		y = numpy.rot90(self.__kartesian)[1]

		r = numpy.rot90(self.__cylinder)[2]
		phi = numpy.rot90(self.__cylinder)[1]

		p_r = numpy.rot90(self.__momentums)[2]
		p_f = numpy.rot90(self.__momentums)[1]

		self.__plot(x, y, 'x(t)', 'y(t)', "$\\vec{r'}$", 'xy_t')
		self.__plot(self.__t, self.__Energy, 't', 'E', 'E(t)', 'E_t')
		self.__plot(self.__t, r, 't', 'r', 'r(t)', 'r_t')
		self.__plot(self.__t, phi, 't', '$\phi$', '$\phi (t)$', 'phi_t')
		self.__plot(self.__t, p_r, 't', '$p_r$', '$p_r (t)$', 'p_r__t')
		self.__plot(self.__t, p_f, 't', '$p_{\phi}$', '$p_{\phi} (t)$', 'p_phi__t')
#----------------------- EXECUTION --------------------------

def Exe(N = N):

	s_1 = [1.5, 1.25*numpy.pi, 0, 0, 1.5**2/2, 0]
	A = charged_particle(N, s_1, 'I', 'darkviolet','o', 'None')
	A.calculate()

	s_2 = [1, 0, 0, 0, -0.5, 0]
	B = charged_particle(N, s_2, 'II', 'red','')
	B.calculate()

	s_3 = [2, 0, 0, 0, -2, 0]
	C = charged_particle(N, s_3, 'III', 'blue','','--', 2)
	C.calculate()

	s_4 = [2, 0, 0, 2, -2, 0]
	D = charged_particle(N, s_4, 'III', 'darkgreen','', 'dotted', 2)
	D.calculate()

	L = [A, B, C, D]

	DRAWING = draw_plots(L)
	DRAWING.plot()

Exe()



