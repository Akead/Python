import sys
import time as TIME
import argparse
import numpy
from scipy.linalg import eigh
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

class Quantum_Oscilator_Galerkin_Method:

    def __init__(self, **kwargs):
        self._L = Atomic_Unit.nano2Bohr(kwargs['L'])
        self._N = kwargs['N']
        assert self._N > 0 and isinstance(self._N, int), f'Invaild N value!'
        self._m = Atomic_Unit.electon_mass()
        self._omega = Atomic_Unit.omega_meV2atomicunits(kwargs['energy'])
    
    def _potential_of_well(self, state):
        assert state <= self._N, f'Invaild state number!'
        return Atomic_Unit.h_bar()**2 * numpy.pi**2 * state**2 / (2 * self._m * self._L**2)
    
    def _eingenfunction_of_potetnial_well(self, i, x):
        return numpy.sqrt(2/self._L) * numpy.sin(i * numpy.pi * x / self._L)

    def _potential(self, i, j):
        assert isinstance(i, int) and i > 0 and isinstance(j, int) and j > 0, f'Error - i or j value is not int!' 
        if i == j:
            return 1/24 * self._L**2 * (i**2 * numpy.pi**2 - 6) * self._m * self._omega**2 / (i**2 * numpy.pi**2)
        else:
            return 2 * self._L**2 * i * j * self._m * self._omega**2 * ((-1)**(i + j) + 1) / (numpy.pi**2 * (-i + j)**2 * (i + j)**2)

    def _hamiltonian_matrix(self):
        kroneker_delta_matrix = numpy.eye(self._N)
        potential = numpy.array([[self._potential(i, j) for j in range(1, self._N + 1)] for i in range(1, self._N + 1)])
        energy =  [self._potential_of_well(i) for i in range(1, self._N + 1)]
        return energy * kroneker_delta_matrix + potential

    def _eingenvalues_and_eingenvectors(self):
        self._eingenvalues, self._eingenvectors = eigh(self._hamiltonian_matrix())

    def return_eingenvalues_in_meV(self, state):
        if not hasattr(self, '_eingenvalues'):
            self._eingenvalues_and_eingenvectors()
        try:
            return Atomic_Unit.Hartree2milieV(self._eingenvalues[state])
        except:
            return None

    def return_eingenvector(self, state):
        if not hasattr(self, 'self._eingenvectors'):
            self._eingenvalues_and_eingenvectors()
        return numpy.rot90(self._eingenvectors)[::-1][state]

    def return_eingenfunction(self, state, N):
        return numpy.array([sum([element * self._eingenfunction_of_potetnial_well(i , x) for i, element in enumerate(self.return_eingenvector(state))]) for x in numpy.linspace(0, self._L, N)])


def zadanie_1(options, L, name, opt = 20):
    print(f'\n{107*"-"} {name} {107*"-"}\n')
    start_time = TIME.time()

    values = []

    for i in range(1,opt + 1):
        configuration = {'L': L, 'N': i, 'energy': 10}
        A = Quantum_Oscilator_Galerkin_Method(**configuration)
        row = [i]
        for state in range(0,6):
            eingenvalue = A.return_eingenvalues_in_meV(state)
            row.append('-' if eingenvalue is None else '{:.7f}'.format(A.return_eingenvalues_in_meV(state)))
        values.append(row)

    table = PrettyTable(['N', *[f'Eingenvalue of H [meV]: State {i}' for i in range(1,7)]])
    for value in values:
        table.add_row([*value])
    table.add_row(['T. value - quantum harmonic oscilator',  *['{:.7f}'.format((i + 0.5) * 10) for i in range(0, 6)]])
    print(table, '\n')

    print(f'{225*"_"}\n')
    return TIME.time() - start_time

def zadanie_2(options, L, name, opt = 20):
    def helper(state):
        colors = ['red', 'coral', 'peru', 'olive', 'lime', 'forestgreen', 'blue', 'darkblue', 'black', 'green', 'indigo']
        eingenfunctions = []
        legends = []
        index = numpy.linspace(state, state + opt, len(colors), dtype = int)
        for i, color in zip(index.tolist(), colors):
            configuration = {'L': L, 'N': 1 + i, 'energy': 10}
            A = Quantum_Oscilator_Galerkin_Method(**configuration)
            eingenfunctions.append(A.return_eingenfunction(state, 200) if A.return_eingenfunction(state, 200)[100] > 0 else -1 * A.return_eingenfunction(state, 200))
            legends.append({'color': color, 'legend': '$\Psi_{}(x)$ -> N: {}'.format(state, i)})

        settings = {'name': 'zadanie_2_{}_nm_N_{}_'.format(L, state), 'title': 'The form of the wave function $\Psi_{}(x)$ for different N'.format(state),
            'x_label': 'x [nm]', 'y_label': '$\Psi_{}(x)$'.format(state), 'location': 'upper right'}
        Plot_1D([*[numpy.linspace(0,100,200) for i in range(len(colors))]], eingenfunctions, [*legends], settings, options['plot'], options['no_save']).plot()

    #print(f'\n{107*"-"} {name} {107*"-"}\n')
    start_time = TIME.time()

    helper(1)
    helper(2)
    helper(3)

    #print(f'{225*"_"}\n')
    return TIME.time() - start_time

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--plot', action='store_true', help='Plot charts show')
    args.add_argument('--no_save', action='store_false', help='Plot charts save')
    args.add_argument('--debug', action='store_true', help='Debuging')
    options = vars(args.parse_args())

    time_table = []
    time_table.append(('1', zadanie_1(options, 100, 'ZADANIE 1')))
    time_table.append(('2', zadanie_2(options, 100, 'ZADANIE 2')))
    time_table.append(('3.1.1', zadanie_1(options, 20, 'ZADANIE 3.1.1', 100)))
    time_table.append(('3.1.2', zadanie_2(options, 20, 'ZADANIE 3.1.2', 100)))
    time_table.append(('3.2.1', zadanie_1(options, 200, 'ZADANIE 3.2.1', 30)))
    time_table.append(('3.2.2', zadanie_2(options, 200, 'ZADANIE 3.2.2', 30)))

    table = PrettyTable(['Task', 'Time [s]'])
    for i, value in time_table:
        table.add_row([i, value])
    print(table, '\n')
