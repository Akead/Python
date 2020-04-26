#!/usr/bin/env python3
import sys
import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



#Chart drawing
draw_chart = (len(sys.argv)<=1)


#-------------------- GLOBAL VARIABLES --------------#

g = 9.81	#m/s^2

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

class bound_particle:

	def __init__(self, alpha, N, t_0, dt, start_cond, name = '', g = g):
		self.__alpha = alpha
		self.__N = N
		self.__name = name
		self.__dt = dt
		self.__g = g
		self.__resoult = [start_cond]
		self.__t = [t_0]
		self.__Energy = numpy.array([])

	def __function(self, s_vector):
		s_0 = s_vector[0]
		s_1 = s_vector[1]
		s_2 = s_vector[2]
		s_3 = s_vector[3]

		d_0 = s_2
		d_1 = s_3
		d_2 = -self.__g * numpy.cos(self.__alpha)**2/numpy.sin(self.__alpha) * numpy.sin(s_0)/s_1 - 2*s_2*s_3/s_1
		d_3 = numpy.sin(self.__alpha)**2*s_1*s_2**2 - self.__g*numpy.sin(self.__alpha)*numpy.cos(self.__alpha)**2 * (1-numpy.cos(s_0))

		return numpy.array([d_0, d_1, d_2, d_3])

	def __cylinder_to_kartesian(self):
		tmp = numpy.rot90(self.__resoult)
		phi = tmp[3]
		z = tmp[2]

		ro = z * numpy.tan(self.__alpha)
		x = ro * numpy.cos(phi)
		y = ro * numpy.sin(phi)
		z = z

		return (x,y,z)

	def __rotate_kartesian(self):
		theta = numpy.pi/2 - self.__alpha
		x, y, z = self.__cylinder_to_kartesian()

		n_x = x*numpy.cos(theta) + z*numpy.sin(theta)
		n_y = y
		n_z = -numpy.sin(theta)*x + z*numpy.cos(theta)

		return (n_x, n_y, n_z)

	def __total_energy(self):
		tmp = numpy.rot90(self.__resoult)
		phi = tmp[3]
		z = tmp[2]
		d_phi = tmp[1]
		d_z = tmp[0]

		self.__Energy = 0.5 * (numpy.tan(self.__alpha)**2 *z**2 *d_phi**2 + (d_z**2)/(numpy.cos(self.__alpha)**2)) + self.__g*z * numpy.sin(self.__alpha) * (1 - numpy.cos(phi))

	def __plot(self,x, y, x_label, y_label, legend_label, name, color = 'blue', draw = draw_chart, save = True):
		fig, ax = plt.subplots()
		ax.plot(x,y)

		#Legend
		legend = []
		description = mlines.Line2D([], [], label = legend_label, color = color)
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

	def __plot_3D(self,x, y, z, x_label, y_label, z_label, legend_label, name, color = 'blue', draw = draw_chart, save = True):
		self.__plot(x, z, "x'", "z'", "z'(x')", name + 'z_x')
		self.__plot(x, y, "x'", "y'", "x'(y')", name + 'x_y')

		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.plot(x,y,z)


		#Legend
		legend = []
		description = mlines.Line2D([], [], label = legend_label, color = color)
		legend.append(description)
		ax.legend(handles = legend, loc = 'lower right', fontsize = 'x-large')

		#Label XYZ
		ax.set(ylabel = y_label, xlabel = x_label, zlabel = z_label)

		#Save
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.zaxis.label.set_size(18)
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

	def print_resoults(self):

		print(80*'-')
		print('phi     ', '	', 'z     ','	', 'd_phi   ', '	','d_z   ', '	', 'E')	
		print(80*'-')
		for i, k in zip(self.__resoult, self.__Energy):
			print('{:.6f}'.format(i[0]), '	','{:.6f}'.format(i[1]), '	', '{:.6f}'.format(i[2]), '	', '{:.6f}'.format(i[3]), '	', '{:.6f}'.format(k), '\n')
		print(80*'-')

	def plot(self):
	
		tmp = numpy.rot90(self.__resoult)
		phi = tmp[3]
		z = tmp[2]
		d_phi = tmp[1]
		d_z = tmp[0]

		self.__plot(self.__t, phi, 't', '$\phi$', '$\phi (t)$', self.__name + 'phi_t')
		self.__plot(self.__t, z, 't', 'z', 'z(t)', self.__name + 'z_t')
		self.__plot(self.__t, d_phi, 't', '$\omega$', '$\omega (t)$', self.__name + 'd_phi_t')
		self.__plot(self.__t, d_z, 't', '$\dot{z}$', '$\dot{z} (t)$', self.__name + 'd_z_t')
		self.__plot(self.__t, self.__Energy, 't', 'E', 'E(t)', self.__name + 'E_t')

		x, y, z = self.__rotate_kartesian()
		self.__plot_3D(x, y, z, "x'", "y'", "z'", "$\\vec{r'}(x',y',z')$", self.__name + '3d')

#----------------------- EXECUTION --------------------------


def Exe(s, a, N = 500):

	A = bound_particle(0.5, N, 0, 0.1, s, a)
	A.calculate()
	A.plot()


s_1 = numpy.array([1.1, 1, 0, 0])
Exe(s_1, 'A')

s_2 = numpy.array([1.1, 0.5, 0, 0])
Exe(s_2, 'B')

s_3 = numpy.array([2.2, 1, 0, 0])
Exe(s_3, 'C')

s_4 = numpy.array([1.1, 1, 0, 5])
Exe(s_4, 'D')

s_5 = numpy.array([1.1, 1, 2, 0])
Exe(s_5, 'F')

