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

class Dispersion:

    def __init__(self, **kwargs):
        self._E = Atomic_Unit.milieV2Hartree(kwargs["E"])
        self._x = Atomic_Unit.nano2Bohr(numpy.array(kwargs['x'], dtype = complex))
        self._delta_x = self._x[1] - self._x[0] 
        self._potential = Atomic_Unit.milieV2Hartree(numpy.array(kwargs['potential'], dtype = complex))
        self._q = self._wave_vector()
        self._wave_function = None

    def _wave_vector(self):
        return numpy.sqrt(2 * Atomic_Unit.electon_mass() * self._E/(Atomic_Unit.h_bar()**2))

    def _wave_function_right(self):
        wave_function = [1+0j, numpy.exp(-1j * self._q * self._delta_x)]
        potential_right = self._potential[::-1]
        for i in range(1, len(potential_right) - 1):
            A = -2 * Atomic_Unit.electon_mass()/(Atomic_Unit.h_bar()**2) * (self._E - potential_right[i]) * self._delta_x**2 * wave_function[i]
            B = -wave_function[i-1] + 2 * wave_function[i]
            wave_function.append(A + B)
        self._wave_function =  numpy.array(wave_function[::-1])

    def wave_function_abs(self):
        self._wave_function_right()
        return abs(self._wave_function)

    def __A(self):
        if self._wave_function is None:
            self._wave_function_right()

        psi_1 = self._wave_function[0]
        psi_2 = self._wave_function[1]
        x_1 = self._x[0]
        x_2 = self._x[1]
        A = psi_1 * numpy.exp(1j * self._q * x_1) - psi_2 * numpy.exp(1j * self._q * x_2)
        B = numpy.exp(1j * self._q * x_1)
        C = numpy.exp(1j * self._q * x_2)
        return  A/(B**2 - C**2)

    def __B(self):
        if self._wave_function is None:
            self._wave_function_right()

        psi_1 = self._wave_function[0]
        psi_2 = self._wave_function[1]
        x_1 = self._x[0]
        x_2 = self._x[1]
        A = (-psi_2 * numpy.exp(1j * self._q * x_1) + psi_1 * numpy.exp(1j * self._q * x_2)) * numpy.exp(1j * self._q * x_1 + 1j * self._q * x_2)
        B = numpy.exp(1j * self._q * x_1)
        C = numpy.exp(1j * self._q * x_2)
        return  -A/(B**2 - C**2)
    
    def return_T(self):
        return 1/abs(self.__A())**2

    def return_R(self):
        return abs(self.__B())**2 / abs(self.__A())**2

    def _wave_function_A_and_B(self, x):
        x = Atomic_Unit.nano2Bohr(x)
        return self.__A() * numpy.exp(1j * self._q * x) + self.__B() * numpy.exp(-1j * self._q * x)

    def wave_function_AB_abs(self, x):
        return abs(self._wave_function_A_and_B(x))

def zadanie_1(options):
    x = numpy.arange(0,340 * 0.5, 0.5)
    potential = numpy.zeros(340)
    potential[100:120] = 20
    potential[220:240] = 20
    configuration = {'x': x, 'potential': potential, 'E': 7}
    Object = Dispersion(**configuration)
    wave_abs = Object.wave_function_abs()
    wave_AB_abs = Object.wave_function_AB_abs(x[:100])
    wave_abs_plot = {'color': 'purple', 'legend': 'Wave function: $\\Psi(x)$'}
    wave_AB_abs_plot = {'color': 'orange', 'legend': 'Wave function: $\\Psi_<(x)$', 'linestyle': 'dashdot'}
    potential_plot = {'color': 'green', 'legend': 'Potential V(x)', 'right_axis': True}

    settings = {'name': 'zadanie_1', 'title': 'Scattering on the double barrier: W = 20 [meV], E = 7 [meV]', 'x_label': 'x [nm]', 'y_label': 'Wave functions: $\\Psi(x)$ and $\\Psi_<(x)$',
    'y_right_label': 'Potential V(x) [meV]', 'location': 'upper right'}
    Plot_1D([x, x[:100], x] , [wave_abs, wave_AB_abs, potential], [wave_abs_plot, wave_AB_abs_plot, potential_plot], settings, options['plot'], options['no_save']).plot()

def zadanie_2(options):
    Energy = numpy.linspace(0, 50, 25000)
    x = numpy.arange(0,340 * 0.5, 0.5)
    potential = numpy.zeros(340)
    potential[100:120] = 20
    potential[220:240] = 20
    T = []
    R = []
    table = PrettyTable(['T(E)', 'E [meV]'])
    for E in Energy:
        configuration = {'x': x, 'potential': potential, 'E': E}
        Object = Dispersion(**configuration)
        t = Object.return_T()
        if t > 0.999:
            table.add_row([t, E])
        T.append(t)
        R.append(Object.return_R())

    print(table, '\n')
    R_plot = {'color': 'purple', 'legend': 'Reflection ratio: $R(E)$'}
    T_plot = {'color': 'orange', 'legend': 'Transmision ratio: $T(E)$', 'right_axis': True}
    settings = {'name': 'zadanie_2', 'title': 'Scattering on the double barrier.\n Reflection $R(E)$ and Transmition $T(E)$ ratio in function of energy E',
    'x_label': 'E [meV]', 'y_label': 'Reflection ratio: $R(E)$', 
    'y_right_label': 'Transmition ratio: $T(E)$', 'location': 'upper right'}
    Plot_1D([Energy, Energy] , [R, T], [R_plot, T_plot], settings, options['plot'], options['no_save']).plot()

    T = [t if t >= 0.9 else None for t in T]

    settings = {'name': 'zadanie_2_1', 'title': 'Scattering on the double barrier.\nTransmition ratio $T(E)$ in function of energy E',
    'x_label': 'E [meV]',
    'y_label': 'Transmition ratio: $T(E)$', 'location': 'upper right'}
    T_plot = {'color': 'orange', 'legend': 'Transmision ratio: $T(E)$'}
    Plot_1D([Energy] , [T], [T_plot], settings, options['plot'], options['no_save']).plot()

def zadanie_3(options):
    Energy = [1.4940597623904956, 5.916236649465979, 13.11852474098964, 22.95691827673107]
    x = numpy.arange(0,340 * 0.5, 0.5)
    potential = numpy.zeros(340)
    potential[100:120] = 20
    potential[220:240] = 20
    for i, E in enumerate(Energy):
        configuration = {'x': x, 'potential': potential, 'E': E}
        Object = Dispersion(**configuration)
        wave_abs = Object.wave_function_abs()
        wave_abs_plot = {'color': 'purple', 'legend': 'Wave function: $\\Psi(x)$'}
        potential_plot = {'color': 'green', 'legend': 'Potential V(x)', 'right_axis': True}
        settings = {'name': f'zadanie_3_{i}', 'title': 'Scattering on the double barrier: W = 20 [meV], E = {:.3f} [meV]'.format(E), 'x_label': 'x [nm]', 'y_label': 'Wave function: $\\Psi(x)$',
        'y_right_label': 'Potential V(x) [meV]', 'location': 'upper right'}
        Plot_1D([x, x] , [wave_abs, potential], [wave_abs_plot, potential_plot], settings, options['plot'], options['no_save']).plot()

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--plot', action='store_true', help='Plot charts show')
    args.add_argument('--no_save', action='store_false', help='Plot charts save')
    args.add_argument('--debug', action='store_true', help='Debuging')
    options = vars(args.parse_args())

    #zadanie_1(options)
    #zadanie_2(options)
    #zadanie_3(options)