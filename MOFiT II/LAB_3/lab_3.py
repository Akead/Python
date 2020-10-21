import sys
import time as TIME
import argparse
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from matplotlib import cm
from copy import deepcopy
from prettytable import PrettyTable

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
    def milieV2joules(unit):
        return unit * 1.602 * 1e-21

    @staticmethod
    def electon_charge():
        return 1

    @staticmethod
    def electon_mass():
        return 0.067
    
    @staticmethod
    def h_bar():
        return 1

    @staticmethod
    def omega_meV2atomicunits(unit):
        return Atomic_Unit.milieV2Hartree(unit)/Atomic_Unit.h_bar()

    @staticmethod
    def omega_atomicunits2meV(unit):
        return Atomic_Unit.Hartree2milieV(unit) * Atomic_Unit.h_bar()

    @staticmethod
    def atomictime2picoseconds(unit):
        return unit * 2.42 * 1e-5

    @staticmethod
    def picoseconds2atomictime(unit):
        return unit/2.42 * 1e5


class Askar_Method_1D:
    '''
    Energy in [meV], metric values in [nm]
    conversion to atomics units. Time in atomic units.
    '''
    def __init__(self, **kwargs):
        self._options = kwargs
        self._omega = Atomic_Unit.omega_meV2atomicunits(kwargs['energy'])
        self._period = 2 * numpy.pi/self._omega
        self._x = Atomic_Unit.nano2Bohr(kwargs['x'])
        self._dx = self._x[1] - self._x[0]
        self._x0 = Atomic_Unit.nano2Bohr(kwargs['x0'])
        self._t = numpy.arange(0, kwargs.setdefault('N_perioids', 10) * self._period, kwargs['dt'])
        self._dt = kwargs['dt']
        self._potential = self._potential_of_harmonic_oscilator() if kwargs.setdefault('potential', None) is None else kwargs['potential']

        self._initialize_eigenfunction_start_values()
        self._eigenfunction_value = []
        self._time = []

    def _potential_of_harmonic_oscilator(self):
        return Atomic_Unit.electon_mass()/2 * self._omega**2 * self._x**2

    def _initialize_eigenfunction_start_values(self):
        self._eigenfunction_value_past = self._eigenfunction_normalize(numpy.exp(-Atomic_Unit.electon_mass()/2 * self._omega * (self._x - self._x0)**2))
        self._eigenfunction_value_past[0] = self._eigenfunction_value_past[-1] = 0+0j
        self._eigenfunction_value_current = self._eigenfunction_normalize(self._eigenfunction_value_past * numpy.exp(-1j/2 * self._omega * self._dt))

    def _eigenfunction_normalize(self, eigenfunction_value):
        I = sum(list(map(lambda i : abs(i)**2 * self._dx, eigenfunction_value)))
        return numpy.array(list(map(lambda i: i/numpy.sqrt(I), eigenfunction_value)))

    def _expected_position_value(self, wave_function, position):
        return Atomic_Unit.Bohr2nano(sum(wave_function.conjugate() * position * wave_function * self._dx).real)

    def _hamiltonian_on_grid(self, wave_function):
        Hamiltonian = [complex(0,0)]
        for i in range(1, len(wave_function)-1):
            A = -Atomic_Unit.h_bar()**2/(2*Atomic_Unit.electon_mass())
            B = wave_function[i+1] + wave_function[i-1] - 2*wave_function[i]
            C = self._potential[i] * wave_function[i]
            Hamiltonian.append(A * B/(self._dx**2)+C)
        Hamiltonian.append(complex(0,0))
        return numpy.array(Hamiltonian)

    def _askar_method(self, wave_function_past, wave_function_curret):
        return wave_function_past + 2 * self._dt/(1j * Atomic_Unit.h_bar()) * self._hamiltonian_on_grid(wave_function_curret)

    def _calculate(self):
        eigenfunction_value_past = self._eigenfunction_value_past
        eigenfunction_value_current = self._eigenfunction_value_current
        self._eigenfunction_value.append(eigenfunction_value_past)
        self._time.append(0.0)
        for i, time in enumerate(self._t):
            eigenfunction_value_future = self._askar_method(eigenfunction_value_past, eigenfunction_value_current)
            eigenfunction_value_past = eigenfunction_value_current
            eigenfunction_value_current = eigenfunction_value_future
            if i and not i%self._options.setdefault('print_step', 100):
                self._eigenfunction_value.append(eigenfunction_value_future)
                self._time.append(time)
            if self._options.setdefault('debug', False) and not i%100:
                print('Calculated {} time steps. This is {:.3f} % of all.'.format(i, 100 * i/len(self._t)))

    def return_probability_density(self):
        if not self._eigenfunction_value:
            self._calculate()
        return (Atomic_Unit.atomictime2picoseconds(numpy.array(self._time)), abs(numpy.array(self._eigenfunction_value)))

    def return_expected_position_value(self):
        if not self._eigenfunction_value:
            self._calculate()
        return (Atomic_Unit.atomictime2picoseconds(numpy.array(self._time)), numpy.array([self._expected_position_value(wave_function, self._x) for wave_function in self._eigenfunction_value]))
        
class Plot_Imshow():

    def __init__(self, data, **kwargs):
        self._data = data[::-1] if kwargs.setdefault('no_y_inverse', True) else data 
        self._kwargs = kwargs

    def _set_figure(self):
        self._figure, self._ax = plt.subplots()
    
    def _show_figure(self):
        self._ax.xaxis.label.set_size(18)
        self._ax.yaxis.label.set_size(18)
        self._ax.title.set_size(18)
        self._ax.tick_params(labelsize = 'large')
        self._ax.get_figure().set_size_inches(18.5, 10.5, forward = True)

        if self._kwargs.setdefault('plot', True):
            plt.show()

    def _save_figure(self, name):

        self._ax.xaxis.label.set_size(18)
        self._ax.yaxis.label.set_size(18)
        self._ax.title.set_size(18)
        self._ax.tick_params(labelsize = 'large')
        self._ax.get_figure().set_size_inches(18.5, 10.5, forward = True)

        if self._kwargs.setdefault('no_save', False):
            self._figure.savefig(f'{name}.png', dpi = 200)
            plt.close()        

    def plot(self):

        self._set_figure()
        imshow = self._ax.imshow(self._data, cmap = self._kwargs.setdefault('heatmap_color', 'jet'), aspect = 'auto')
        
        if not self._kwargs.setdefault('no_colorbar', False):
            colorbar = self._figure.colorbar(imshow)
            colorbar.ax.tick_params(labelsize = 18)
            if self._kwargs.setdefault('colorbar_title', None) is not None:
                colorbar.set_label(self._kwargs['colorbar_title'], size=18)
 
        if self._kwargs.setdefault('extent', None) is not None:
            imshow.set_extent(self._kwargs['extent'])

        self._ax.set_title(self._kwargs.setdefault('title', ''))
        self._ax.set_xlabel(self._kwargs.setdefault('xlabel', ''))
        self._ax.set_ylabel(self._kwargs.setdefault('ylabel', ''))
        self._ax.grid()
        self._show_figure()
        self._save_figure(self._kwargs['name'])

class Plot_1D:

    def __init__(self, x, y, plot_options, settings, plot = True, save = True):
        self._x = x if isinstance(x, list) else [x]
        self._y = y if isinstance(x, list) else [y]
        self._plot_options = plot_options if isinstance(plot_options, list) else [y]
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



def zadanie_1(options):
    print(f'\n{50*"-"} ZADANIE 1 {50*"-"}\n')
    start_time = TIME.time()
    configuration = {'x': numpy.linspace(-100,100,101, dtype = complex), 'energy': 10, 'dt': 2, 'x0':20, 'debug': options['debug']}
    Askar = Askar_Method_1D(**configuration)
    time, probability_density = Askar.return_probability_density()
    options = {'name' : 'zadanie_1', 'title' : 'Time and position dependence of the probability density function $|\psi(x,y)|^2$','xlabel': 'Position x [nm]', 'ylabel': 'Time t [ps]' ,
            'no_save': options['no_save'], 'plot': options['plot'],
            'colorbar_title': '$|\psi(x,y)|^2$', 'extent': [-100, 100, 0, max(time)], 'heatmap_color': 'plasma'}
    Plot_Imshow(probability_density, **options).plot()
    print(f'{111*"_"}\n')
    return TIME.time() - start_time

def zadanie_2(options):

    print(f'\n{50*"-"} ZADANIE 2 {50*"-"}\n')
    start_time = TIME.time()
    configuration = {'x': numpy.linspace(-100,100,101, dtype = complex), 'energy': 10, 'dt': 2, 'x0':20, 'debug': options['debug']}
    Askar = Askar_Method_1D(**configuration)
    time, expected_position_value = Askar.return_expected_position_value()
    x = Atomic_Unit.Bohr2nano(Askar._x0) * numpy.cos(Askar._omega * Atomic_Unit.picoseconds2atomictime(time))
    quantum = {'color': 'purple', 'legend': 'Quantum oscilator: <x>(t)', 'marker': 'o'}
    classic = {'color': 'orange', 'legend': 'Classic oscilator: x(t)', 'marker': '2'}

    settings = {'name': 'zadanie_2', 'title': 'Dependence of the expected position \n value <x>(t) on time for a quantum oscillator and x(t) for a classical oscillator',
            'x_label': 'Time t [ps]', 'y_label': 'Position <x>(t)', 'location': 'upper right'}

    Plot_1D([time, time], [expected_position_value, x], [quantum, classic], settings, options['plot'], options['no_save']).plot()
    print(f'{111*"_"}\n')
    return TIME.time() - start_time


def zadanie_3(options):
    print(f'\n{50*"-"} ZADANIE 3 {50*"-"}\n')
    start_time = TIME.time()
    configuration = {'x': numpy.linspace(-100,100,101, dtype = complex), 'energy': 10, 'dt': 2, 'x0':0, 'debug': options['debug']}
    Askar = Askar_Method_1D(**configuration)
    time, probability_density = Askar.return_probability_density()
    options = {'name' : 'zadanie_3', 'title' : 'Time and position dependence of the probability density function $|\psi(x,y)|^2$','xlabel': 'Position x [nm]', 'ylabel': 'Time t [ps]' ,
            'no_save': options['no_save'], 'plot': options['plot'],
            'colorbar_title': '$|\psi(x,y)|^2$', 'extent': [-100, 100, 0, max(time)], 'heatmap_color': 'plasma'}
    Plot_Imshow(probability_density, **options).plot()
    print(f'{111*"_"}\n')
    return TIME.time() - start_time

def zadanie_4(options):
    print(f'\n{50*"-"} ZADANIE 4 {50*"-"}\n')
    start_time = TIME.time()
    configuration = {'x': numpy.linspace(-100,100,101, dtype = complex), 'energy': 10, 'potential': numpy.zeros(101),'dt': 2, 'x0':0, 'debug': options['debug']}
    Askar = Askar_Method_1D(**configuration)
    time, probability_density = Askar.return_probability_density()
    options = {'name' : 'zadanie_4', 'title' : 'Time and position dependence of the probability density function $|\psi(x,y)|^2$','xlabel': 'Position x [nm]', 'ylabel': 'Time t [ps]' ,
            'no_save': options['no_save'], 'plot': options['plot'],
            'colorbar_title': '$|\psi(x,y)|^2$', 'extent': [-100, 100, 0, max(time)], 'heatmap_color': 'plasma'}
    Plot_Imshow(probability_density, **options).plot()
    print(f'{111*"_"}\n')
    return TIME.time() - start_time

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--plot', action='store_true', help='Plot charts show')
    args.add_argument('--no_save', action='store_false', help='Plot charts save')
    args.add_argument('--debug', action='store_true', help='Debuging')
    options = vars(args.parse_args())

    time_table = []
    time_table.append((1, zadanie_1(options)))
    time_table.append((2, zadanie_2(options)))
    time_table.append((3, zadanie_3(options)))
    time_table.append((4, zadanie_4(options)))

    table = PrettyTable(['Task', 'Time [s]'])
    for i, value in time_table:
        table.add_row([i, value])
    print(table, '\n')
