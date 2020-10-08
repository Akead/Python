import sys
import time
import argparse
import random
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from copy import deepcopy

class Atomic_Unit:
	'''
	Conversion unit class
	'''
	@staticmethod
	def nano2Bohr(unit):
		return unit/5.2917721092 * 1e2

	@staticmethod
	def Bohr2nano(unit):
		return unit * 5.2917721092 * 1e-2

	@staticmethod
	def milieV2Hartree(unit):
		return unit/27.2114 * 1e-3

	@staticmethod
	def Hartree2milieV(unit):
		return unit * 27.2114 * 1e3

	@staticmethod
	def electon_charge():
		return 1

	@staticmethod
	def electon_mass():
		return 0.067
	
	@staticmethod
	def h_bar():
		return 1

class Shot_Method:

	'''
	Energy in [meV], metric values in [nm]
	conversion to atomics units.
	'''

	def __init__(self, configuration, plot_options = None):
		self._options = deepcopy(configuration)
		self._options['dx_nm'] = configuration['L']/configuration['N']
		self._options['dx'] = Atomic_Unit.nano2Bohr(self._options['dx_nm'])
		self._options['E'] = Atomic_Unit.milieV2Hartree(configuration['E'])
		self._potential = Atomic_Unit.milieV2Hartree(numpy.array(self._options['potential']))
		self._results = {}
		self._results[plot_options] = plot_options
		self._initialize_eigenfunction_value()

	def _eigenfunction_normalize(self, eigenfunction_value):
		I = sum(list(map(lambda i : abs(i)**2 * self._options['dx'], eigenfunction_value)))
		return list(map(lambda i: i/numpy.sqrt(I), eigenfunction_value))

	def _initialize_eigenfunction_value(self):
		self._eigenfunction_value = numpy.zeros(self._options['N'])
		self._eigenfunction_value[0] = 0
		self._eigenfunction_value[1] = 1

	def _eigenfunction_value_next(self, value_back, value_current, potential = 0):
			return -2 * Atomic_Unit.electon_mass()/(Atomic_Unit.h_bar()**2) * (self._options['E'] - potential) * self._options['dx']**2 * value_current - value_back + 2*value_current

	def _calculate(self):
		for i in range(1, self._options['N']-1):
			value_back = self._eigenfunction_value[i-1]
			value_current = self._eigenfunction_value[i]
			potential = self._potential[i]
			self._eigenfunction_value[i+1] = (self._eigenfunction_value_next(value_back, value_current, potential))
		self._results['wave_function'] = self._eigenfunction_normalize(self._eigenfunction_value)

	def return_wave_function(self):
		self._calculate()
		return self._results


class Plot_1D:

	def __init__(self, x, y, plot_options, settings, plot = True, save = True):
		self._x = x if isinstance(x, list) else [x]
		self._y = y if isinstance(x, list) else [y]
		self._plot_options = plot_options if isinstance(plot_options, list) else [plot_options]
		self._settings = settings
		self._settings['name'] = f'{settings["name"]}.png'
		self._plot = plot
		self._save = save

	def plot(self):
		
		fig, ax = plt.subplots()
		legend = []
		right_flag = False
		for x, y, options in zip(self._x, self._y, self._plot_options):
			if self._settings.setdefault('log', False):
				if options.setdefault('right_axis', None) is not None:
					right_flag = True
					ax_right = ax.twinx()
					line, = ax_right.loglog(x,y)
				else:
					right_flag = False
					line, = ax.loglog(x, y)

			else:
				if options.setdefault('right_axis', None) is not None:
					right_flag = True
					ax_right = ax.twinx()
					line, = ax_right.plot(x,y)
				else:
					right_flag = False
					line, = ax.plot(x, y)

			if self._settings.setdefault('xticks', None) is not None:
				ax.set_xticks(self._settings['xticks'])
				if self._settings.setdefault('xticks_labels', None) is not None:
					ax.set_xticklabels(self._settings['xticks_labels'])

			color = options['color'] if options.setdefault('color', None) is not None else 'black'
			linestyle = options['linestyle'] if options.setdefault('linestyle', None) is not None else 'solid'
			marker = options['marker'] if options.setdefault('marker', None) is not None else None
			markersize = options['markersize'] if options.setdefault('markersize', None) is not None else 6
			label = options['legend'] if options.setdefault('legend', None) is not None else ''
			line.set_color(color)
			line.set_linestyle(linestyle)
			line.set_marker(marker)
			line.set_markersize(markersize)
			description = mlines.Line2D([], [], label = label, color = color, marker = marker, linestyle = linestyle)
			legend.append(description)

		location = self._settings['location'] if self._settings.setdefault('location', None) is not None else 'lower right'
		if label:
			ax.legend(handles = legend, loc = location, fontsize = 'large')
		ax.set(title = self._settings['title'], xlabel = self._settings['x_label'], ylabel = self._settings['y_label'])
		if right_flag:
			ax_right.set(ylabel = self._settings['y_right_label'])
			ax_right.yaxis.label.set_size(18)
			ax_right.tick_params(labelsize = 'large')
			ax_right.grid()
		ax.grid()
		ax.xaxis.label.set_size(18)
		ax.yaxis.label.set_size(18)
		ax.title.set_size(18)
		ax.tick_params(labelsize = 'large')
		ax.get_figure().set_size_inches(18.5, 10.5, forward = True)

		if self._plot:
			plt.show()
		
		if self._save:
			fig.savefig(self._settings['name'], dpi = 200)
			plt.close()

def zero_places_EV(n, L=100):
	return Atomic_Unit.Hartree2milieV(n**2  * Atomic_Unit.h_bar()**2 * numpy.pi**2 / (2 * Atomic_Unit.electon_mass() * Atomic_Unit.nano2Bohr(L)**2))

def zad_1():
	print(f'\n{30*"-"} Zadanie 1 {30*"-"}\n')
	potential = numpy.zeros(101)
	energy = numpy.arange(0, 35, 0.1)
	function_value = []
	for E in energy:
		configuration = {'N': 101, 'L': 100,  'E': E, 'potential': potential}
		Object = Shot_Method(configuration)
		resoults =  Object.return_wave_function()
		wave_function = resoults['wave_function']
		function_value.append(wave_function[-1])
		if options['debug']:
			print(f'Calculated for energy E: {len(function_value)/10 - 1} meV')
			print(f'Function value: {function_value[int(10*E)]}')
	
	zero_places_y = numpy.zeros(7)
	zero_places_x = list(map(zero_places_EV, range(1,8)))

	wave_function_plot = {'color': 'red', 'legend': '$\Psi_N(E)$'}
	zero_line_plot = {'color': 'blue', 'legend': 'y(E) = 0 '}
	zero_places_plot = {'color': 'green', 'marker': 'd', 'markersize': 10, 'linestyle': 'None', 'legend': 'Wartość dokładna: $E_n$'}
	settings = {'name': 'zad_1_1', 'title': 'Wykres zależności wartości funkcji falowej $\Psi_N(E)$ na prawym brzegu studni potencjału od energii $E$', 'x_label': 'E [meV]', 'y_label': '$\Psi_N(E)$', 'location': 'upper right'}
	Plot_1D([energy, [0,35], zero_places_x], [function_value, [0,0], zero_places_y], [wave_function_plot, zero_line_plot, zero_places_plot], settings, options['plot'], options['no_save']).plot()

	selected_energy = [27.9, 29.3, 26.5]
	file_names = ['2_1', '2_2', '2_3']

	for E, label in zip(selected_energy, file_names):
		potential = numpy.zeros(501)
		configuration = {'N': 501, 'L': 100,  'E': E, 'potential': potential}
		Object = Shot_Method(configuration)
		resoults =  Object.return_wave_function()
		wave_function = resoults['wave_function']
		wave_function_x = numpy.arange(501)
		wave_function_plot = {'color': 'green'}
		settings = {'name': f'zad_1_{label}', 'title': f'Wykres funkcji falowej $\Psi(x)$ dla energii $E = {E}$ [meV]', 'x_label': 'x [nm]', 'y_label': '$\Psi(x)$', 'xticks': numpy.arange(0, 501,50), 'xticks_labels': numpy.arange(0, 101,10)}
		Plot_1D([wave_function_x], [wave_function], [wave_function_plot], settings, options['plot'], options['no_save']).plot()

def list_bisection(value_a, value_b, np_array_x, np_array_y):

	a = numpy.where(numpy.isclose(np_array_x, value_a, 1e-6))[0][0]
	b = numpy.where(numpy.isclose(np_array_x, value_b, 1e-6))[0][0]


	while True:

		x = (a+b)//2
		X = np_array_y[x]
		A = np_array_y[a]
		B = np_array_y[b]

		if abs(a-b) <=1 :
			return(x, X)

		elif X * B < 0:
			a = x

		elif X * A < 0:
			b = x

		else:
			print('Bisection method error')
			assert 0, f'Bisection method error'



def zad_2():
	def helper(N, dE):
		potential = numpy.zeros(N)
		function_value = []
		for E in energy:
			configuration = {'N': N, 'L': 100,  'E': E, 'potential': potential}
			Object = Shot_Method(configuration)
			resoults =  Object.return_wave_function()
			wave_function = resoults['wave_function']
			function_value.append(wave_function[-1])
			if options['debug']:
				print(f'Calculated for energy E: {len(function_value)*dE -dE} meV')
				print(f'Function value: {function_value[int(E/dE)]}')

		range_list = [(0.0, 1.0), (2.0, 3.0), (4.0, 6.0), (7.0, 10.0), (13, 16), (19, 21), (26, 29)]
		bisection_values_energy = []
		bisection_values_function = []

		for scope in range_list:
			E, function = list_bisection(*scope, energy, function_value)
			bisection_values_energy.append(E * dE)
			bisection_values_function.append(function)

		print(f'\n\n{30*"_"} Zadanie 2: N = {N} {30*"_"}\n')
		print(f'\n{80*"_"}\n')
		print(f'Bisection method: E [meV] | Bisection method: F | Extract value: E [meV] | State')
		print(f'{80*"-"}')
		for bis_E, bis_F, ext_E, number in zip(bisection_values_energy, bisection_values_function, zero_places_x, range(len(zero_places_x))):
			print('{:.3f} \t\t\t  | {:f} \t\t| {:f} \t\t | {}'.format(bis_E, bis_F, ext_E, number))
		print(f'{80*"-"}\n')

		

	print(f'\n{30*"-"} Zadanie 2 {30*"-"}\n')
	zero_places_x = list(map(zero_places_EV, range(0,8)))

	dE = 0.001
	energy = numpy.arange(0, 35, dE)
	helper(101, dE)

	dE = 0.001
	energy = numpy.arange(0, 35, dE)
	helper(301, dE)

def zad_3():

	def helper(N, dE, W, range_list, i = 0):
		potential[N//2] = -W
		function_value = []
		for E in energy:
			configuration = {'N': N, 'L': 100,  'E': E, 'potential': potential}
			Object = Shot_Method(configuration)
			resoults =  Object.return_wave_function()
			wave_function = resoults['wave_function']
			function_value.append(wave_function[-1])

		bisection_values_energy = []
		bisection_values_function = []


		for scope in range_list:
			E, function = list_bisection(*scope, energy, function_value)
			bisection_values_energy.append(E * dE -50 )
			bisection_values_function.append(function)

		if options['debug']:
			print(f'\n{80*"_"}\n')
			print(f'Bisection method: E [meV] | Bisection method: F | State')
			print(f'{80*"-"}')
			for bis_E, bis_F, number in zip(bisection_values_energy, bisection_values_function,  range(i,8)):
				print('{:.3f} \t\t\t  | {:f} \t\t| {}'.format(bis_E, bis_F, number))
			print(f'{80*"-"}\n')

		return bisection_values_energy

	print(f'\n{30*"-"} Zadanie 3 {30*"-"}\n')
	N = 301
	dE = 0.1
	energy = numpy.arange(-50, 50, dE)
	potential = numpy.zeros(N)
	potential_level_zero = numpy.arange(0, 1001, 100)

	zero_places_zero = []

	range_list = [(-50, 0.6)]

	for W in potential_level_zero:
		print(f'\n\n{30*"_"} Zadanie 3: W = {W} {30*"_"}\n')

		zero_places_zero.append(*helper(N, dE, W, range_list))

	potential_level = numpy.arange(0, 3001, 100)
	zero_places = []
	range_list = [(2.0, 2.3), (2.3, 6.0), (8.0, 9.1), (9.1, 20.0), (20.0, 21.0), (21.0, 35.0), (30.0, 37.0)]

	for W in potential_level:
		print(f'\n\n{30*"_"} Zadanie 3: W = {W} {30*"_"}\n')

		zero_places.append(helper(N, dE, W, range_list, 1))
	zero_places = numpy.rot90(numpy.array(zero_places))[::-1]
	zero_places_x = numpy.array([potential_level, potential_level, potential_level, potential_level, potential_level, potential_level, potential_level])
	plot_options = [
		{'color': 'red', 'legend': 'Miejsce zerowe stanu: 1', 'marker': 'd'},
		{'color': 'blue', 'legend': 'Miejsce zerowe stanu: 2', 'marker': 'o'},
		{'color': 'green', 'legend': 'Miejsce zerowe stanu: 3', 'marker': '*'},
		{'color': 'coral', 'legend': 'Miejsce zerowe stanu: 4', 'marker': 'd'},
		{'color': 'peru', 'legend': 'Miejsce zerowe stanu: 5', 'marker': 'o'},
		{'color': 'forestgreen', 'legend': 'Miejsce zerowe stanu: 6', 'marker': '*'},
		{'color': 'darkblue', 'legend': 'Miejsce zerowe stanu: 7', 'marker': 'd'},
		{'color': 'black', 'legend': 'Miejsce zerowe stanu: 0', 'marker': 'o'}
		]

	settings = {'name': 'zad_3_1', 'title': 'Wykres zależności miejsc zerowych w funkcji W [meV]', 'x_label': 'W [meV]', 'y_label': 'E [meV]', 'location': 'lower right'}
	Plot_1D([*zero_places_x, potential_level_zero], [*zero_places, zero_places_zero], [*plot_options], settings, options['plot'], options['no_save']).plot()



def zad_3_1():

	W = 1000
	potential = numpy.zeros(301)
	potential[301//2] = -W
	energy = [2.250, 2.590, 9.030, 10.320, 20.330, -48.4152, ]
	x = numpy.arange(0,100,100/301)

	function_value = []
	for E in energy:
		configuration = {'N': 301, 'L': 100,  'E': E, 'potential': potential}
		Object = Shot_Method(configuration)
		resoults =  Object.return_wave_function()
		wave_function = resoults['wave_function']
		function_value.append(wave_function)
		if options['debug']:
			print(f'Calculated for energy E: {E} meV')
			print(f'Function value: {function_value[-1][-1]}')

	plot_options = [
		{'color': 'red', 'legend': 'Stan: 1'},
		{'color': 'blue', 'legend': 'Stan: 2'},
		{'color': 'green', 'legend': 'Stan: 3'},
		{'color': 'coral', 'legend': 'Stan: 4'},
		{'color': 'peru', 'legend': 'Stan: 5'},
		{'color': 'black', 'legend': 'Stan podstawowy: 0'}
		]

	settings = {'name': 'zad_3_2_1', 'title': 'Wykres poszczególnych stanów funkcji falowej $\Psi(x)$', 'x_label': 'x [nm]', 'y_label': '$\Psi(x)$', 'location': 'upper left'}
	Plot_1D([x,x,x,x,x,x], [*function_value], [*plot_options], settings, options['plot'], options['no_save']).plot()

	settings = {'name': 'zad_3_2_2', 'title': 'Wykres poszczególnych stanów funkcji falowej $\Psi(x)$', 'x_label': 'x [nm]', 'y_label': '$\Psi(x)$', 'location': 'upper left'}
	Plot_1D([x,x,x,x,x], [*function_value[:-1]], [*plot_options[:-1]], settings, options['plot'], options['no_save']).plot()

def zad_3_2():

	def helper(energy, colors, file_name, title):
		potential_level = numpy.arange(0, 3001, 300)
		potential = numpy.zeros(301)
		x = numpy.arange(0,100,100/301)

		function_value = []

		for E, W in zip(energy, potential_level):
			potential[301//2] = -W
			configuration = {'N': 301, 'L': 100,  'E': E, 'potential': potential}
			Object = Shot_Method(configuration)
			resoults =  Object.return_wave_function()
			wave_function = resoults['wave_function']
			function_value.append(wave_function)
			if options['debug']:
				print(f'Calculated for energy E: {E} meV')

		X = [x for _ in range(len(function_value))]

		plot_options = [{'color': color, 'legend': 'Funkcja falowa $\Psi(x)$ dla W: {:d} [meV]'.format(label)} for color, label in zip(colors, potential_level)]
		settings = {'name': f'zad_3_3_{file_name}', 'title': f'{title}', 'x_label': 'x [nm]', 'y_label': '$\Psi(x)$', 'location': 'lower left'}
		Plot_1D([*X], [*function_value], [*plot_options], settings, options['plot'], options['no_save']).plot()

	energy_1 = [2.250 for _ in range(0,31)]
	energy_2 = [5.080, 4.410, 3.830, 3.420, 3.140, 2.960, 2.840, 2.750, 2.690, 2.630, 2.590, 2.560, 2.530, 2.510, 2.490, 2.480, 2.460, 2.450, 2.440, 2.430, 2.420, 2.410, 2.400, 
				2.400, 2.390, 2.380, 2.380, 2.370, 2.370, 2.370, 2.360]
	
	E_1 = energy_1[::3]
	E_2 = energy_2[::3]
	colors = ['blueviolet', 'crimson', 'purple', 'navy', 'dodgerblue', 'lime', 'yellow', 'red', 'green', 'orange', 'black']

	helper(E_1, colors, 1, 'Wykres zależności funkcji falowej $\Psi(x)$ od W [meV] dla stanu: 1')
	helper(E_2, colors, 2, 'Wykres zależności funkcji falowej $\Psi(x)$ od W [meV] dla stanu: 2')



if __name__ == "__main__":
	args = argparse.ArgumentParser()
	args.add_argument('--plot', action='store_true', help='Plot charts show')
	args.add_argument('--no_save', action='store_false', help='Plot charts save')
	args.add_argument('--debug', action='store_true', help='Debuging')
	options = vars(args.parse_args())

	#zad_1()
	#zad_2()
	#zad_3()
	#zad_3_1()
	#zad_3_2()