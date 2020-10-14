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
    def electon_charge():
        return 1

    @staticmethod
    def electon_mass():
        return 0.067
    
    @staticmethod
    def h_bar():
        return 1

class Imagination_Time_Method_1D:
    '''
    Energy in [meV], metric values in [nm]
    conversion to atomics units.
    '''
    def __init__(self, alpha, configuration, plot_options, force_iter = False):
        self._options = deepcopy(configuration)
        self._options['dx'] = Atomic_Unit.nano2Bohr(self._options['dx'])
        self._options['dx_nm'] = configuration['dx']
        self._potential = Atomic_Unit.milieV2Hartree(numpy.array(self._options['potential'], dtype=complex))
        self._alpha = alpha*0.95*Atomic_Unit.electon_mass()*self._options['dx']**2/Atomic_Unit.h_bar()**2
        self._alpha_val = alpha
        self._iteration = configuration.setdefault('iteration', None)
        self._options['precision'] = self._options['precision'] if self._iteration is None else -1
        self._results = {}
        self._plot_options = plot_options
        self._force_iter = force_iter

    def _initialize_eigenfunction_value(self):
        eigenfunction_value = numpy.ones(self._options['N'], dtype=complex)
        eigenfunction_value = list(map(lambda i: complex(random.uniform(-1,1), random.uniform(-1,1)) , eigenfunction_value))
        eigenfunction_value[0] = eigenfunction_value[-1] = complex(0,0)
        return self._eigenfunction_normalize(eigenfunction_value)

    def _eigenfunction_normalize(self, eigenfunction_value):
        I = sum(list(map(lambda i : abs(i)**2 * self._options['dx'], eigenfunction_value)))
        return list(map(lambda i: i/numpy.sqrt(I), eigenfunction_value))

    def _expected_energy_value(self, wave_function, hamiltonian):
        return sum(list(map(lambda x, y: x.conjugate() * y* self._options['dx'], wave_function, hamiltonian))).real

    def _hamiltonian_on_grid(self, wave_function):
        Hamiltonian = [complex(0,0)]
        for i in range(1, len(wave_function)-1):
            A = -Atomic_Unit.h_bar()**2/(2*Atomic_Unit.electon_mass())
            B = wave_function[i+1] + wave_function[i-1] - 2*wave_function[i]
            C = self._potential[i] * wave_function[i]
            Hamiltonian.append(A * B/(self._options['dx']**2)+C)
        Hamiltonian.append(complex(0,0))
        return Hamiltonian

    def _basic_state_calculation(self):
        self._results['basic_state'] = {'energy': [], 'wave_function': None, 'plot': self._plot_options['basic_state']}
        self._results['basic_state']['wave_function'] = self._initialize_eigenfunction_value()
        while len(self._results['basic_state']['energy']) <= 1 or abs((self._results['basic_state']['energy'][-2] - self._results['basic_state']['energy'][-1]))/(self._results['basic_state']['energy'][-2]) > self._options['precision'] or self._force_iter:
            hamiltonian = self._hamiltonian_on_grid(self._results['basic_state']['wave_function'])
            self._results['basic_state']['energy'].append(self._expected_energy_value(self._results['basic_state']['wave_function'], hamiltonian))
            self._results['basic_state']['wave_function'] = self._eigenfunction_normalize(list(map(lambda x, y: x - self._alpha*y, self._results['basic_state']['wave_function'], hamiltonian)))
            if self._iteration is not None and self._iteration < len(self._results['basic_state']['energy']):
                break
            if options.setdefault('debug', False) and not len(self._results['basic_state']['energy'])%50:
                print(f"Basic state iteration: {len(self._results['basic_state']['energy'])}")	

    def _first_state_calculation(self):
        if self._results.setdefault('basic_state', None) is None:
            print(f"Basic state not defined.")
            self._basic_state_calculation()		
        self._results['first_state'] = {'energy': [], 'wave_function': None, 'plot': self._plot_options['first_state']}
        self._results['first_state']['wave_function'] = self._initialize_eigenfunction_value()
        while len(self._results['first_state']['energy']) <= 1 or abs((self._results['first_state']['energy'][-2] - self._results['first_state']['energy'][-1]))/(self._results['first_state']['energy'][-2]) > self._options['precision'] or self._force_iter:
            hamiltonian = self._hamiltonian_on_grid(self._results['first_state']['wave_function'])
            self._results['first_state']['energy'].append(self._expected_energy_value(self._results['first_state']['wave_function'], hamiltonian))
            self._results['first_state']['wave_function'] = self._eigenfunction_normalize(list(map(lambda x, y: x - self._alpha*y, self._results['first_state']['wave_function'], hamiltonian)))
            C = sum(list(map(lambda x, y: x.conjugate() * y * self._options['dx'], self._results['basic_state']['wave_function'], self._results['first_state']['wave_function'])))
            if options.setdefault('debug', False) and not len(self._results['first_state']['energy'])%50:
                print(f"First state iteration: {len(self._results['first_state']['energy'])}")
            self._results['first_state']['wave_function'] = self._eigenfunction_normalize(list(map(lambda x, y: x - C * y, self._results['first_state']['wave_function'], self._results['basic_state']['wave_function'])))
            if self._iteration is not None and self._iteration < len(self._results['first_state']['energy']):
                break

    def _x_axis(self):
        return [self._options['dx_nm'] * i for i in range(self._options['N'])]

    def basic_state(self):
        if self._results.setdefault('basic_state', None) is None:
            self._basic_state_calculation()
        return self.__output('basic_state')

    def first_state(self):
        if self._results.setdefault('basic_state', None) is None:
            self._basic_state_calculation()
        if self._results.setdefault('first_state', None) is None:
            self._first_state_calculation()
        return self.__output('first_state')

    def __output(self, state):
        self._plot_options[state]['legend'] = None
        list_N = [i for i in range(len(self._results[state]['energy']))]
        self._plot_options[state]['legend'] = '{} state: $\\alpha = {:.2f} = {:.2f} \\alpha_o$ \n Iteration: {} '.format(state[:-6].capitalize(), self._alpha, self._alpha_val, len(list_N))
        return {'wave_function': {'x': self._x_axis() ,'y': self._results[state]['wave_function']},
                'energy': {'x': list_N, 'y': self._results[state]['energy']}, 'plot_options': self._plot_options[state]}

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

def zad_1(parameters, colors, options):
    configuration = {'N': 300, 'dx': 100/300, 'precision': 1e-6, 'debug': options['debug']}
    configuration['potential'] = numpy.zeros(configuration['N'])
    Objects = [Imagination_Time_Method_1D(parameter, configuration,
            {'basic_state': {'color': color['basic_state']['color']},
             'first_state': {'color': color['first_state']['color']}})
              for parameter, color in zip(parameters, colors)]

    X = []
    Y = []
    plot_options = []
    for obj in Objects:
        results = obj.basic_state()
        X.append(numpy.array(results['energy']['x']))
        Y.append(Atomic_Unit.Hartree2milieV(numpy.array(results['energy']['y'])))
        plot_options.append(results['plot_options'])
    settings = {'name': 'zad_1_1e-6', 'title': 'Convergence of basic state <E> depending on value of $\\alpha$', 'x_label': 'Iteration', 'y_label': 'Basic state <E>', 'log': True, 'location': 'upper right'}
    Plot_1D(X, Y, plot_options, settings, options['plot'], options['no_save']).plot()

    table = PrettyTable(['State', 'Alpha value', '<E> [meV]', 'Iterations'])

    energy = [i[-1] for i in Y]
    iterations = [i[-1] for i in X]
    for alpha, E, i in zip(parameters, energy, iterations):
        table.add_row(['Basic','{:.2f}'.format(alpha), E, i])
    print(table, '\n')

def zad_2(options):
    configuration = {'N': 300, 'dx': 100/300, 'precision': 1e-7, 'debug': options['debug']}
    configuration['potential'] = numpy.zeros(configuration['N'])
    Object = Imagination_Time_Method_1D(0.95, configuration, {'basic_state': None, 'first_state': {'color': 'darkblue'}})
    results = Object.first_state()
    wave_function_x = numpy.array(results['wave_function']['x'])
    wave_function_y = numpy.array(results['wave_function']['y'])
    wave_function_plot = results['plot_options']
    wave_function_plot['legend'] = 'First state wave function'
    potential_y = numpy.array(configuration['potential']).real
    potential_x = numpy.linspace(0,100,300)
    potential_plot = {'color': 'green', 'legend': 'Potential V(x)', 'right_axis': True}
    wave_function_plot_abs = {'color': 'orange', 'legend': 'First state wave function module'}
    settings = {'name': 'zad_2_wave', 'title': 'Wave function of first state', 	'x_label': 'x [nm]', 'y_label': 'First state wave function: $\\Psi_1(x)$', 'y_right_label': 'Potential V(x) [meV]', 'location': 'lower right'}
    Plot_1D([wave_function_x.real, wave_function_x.real, potential_x] , [wave_function_y.real, abs(wave_function_y), potential_y], [wave_function_plot, wave_function_plot_abs, potential_plot], settings, options['plot'], options['no_save']).plot()

    energy_x = numpy.array(results['energy']['x']).real
    energy_y = Atomic_Unit.Hartree2milieV(numpy.array(results['energy']['y']).real)
    energy_plot = results['plot_options']
    energy_plot['legend'] = None
    settings = {'name': 'zad_2_convergence', 'title': 'Documentation of the convergence of the solution for first state','x_label': 'Iteration', 'y_label': 'First state <E>', 'log': True}
    Plot_1D(energy_x, energy_y, [energy_plot], settings, options['plot'], options['no_save']).plot()

    table = PrettyTable(['State', '<E> [meV]'])
    table.add_row(['First', energy_y[-1]])
    print(table, '\n')

def zad_3(options):
    configuration = {'N': 300, 'dx': 100/300, 'iteration': 2e4, 'precision': 0, 'debug': options['debug']}
    configuration['potential'] = numpy.zeros(301)
    configuration['potential'][150] = -2000
    Object = Imagination_Time_Method_1D(0.95, configuration, {'basic_state': {'color': 'red'}, 'first_state': {'color': 'darkblue'}}, force_iter=True)
    results = Object.first_state()
    wave_function_first_x = numpy.array(results['wave_function']['x'])
    wave_function_first_y = numpy.array(results['wave_function']['y'])
    wave_function_plot_first = results['plot_options']
    wave_function_plot_first['legend'] = 'First state wave function'
    energy_x_first = numpy.array(results['energy']['x']).real
    energy_y_first = Atomic_Unit.Hartree2milieV(numpy.array(results['energy']['y']).real)
    energy_plot_first = results['plot_options']
    energy_plot_first['legend'] = 'First state'

    results = Object.basic_state()
    wave_function_basic_x = numpy.array(results['wave_function']['x'])
    wave_function_basic_y = numpy.array(results['wave_function']['y'])
    wave_function_plot_basic = results['plot_options']
    wave_function_plot_basic['legend'] = 'Basic state wave function'
    energy_x_basic = numpy.array(results['energy']['x']).real
    energy_y_basic = Atomic_Unit.Hartree2milieV(numpy.array(results['energy']['y']).real)
    energy_plot_basic = results['plot_options']
    energy_plot_basic['legend'] = 'Basic state'

    potential_y = numpy.array(configuration['potential']).real
    potential_x = numpy.linspace(0,100,301)
    potential_plot = {'color': 'green', 'legend': 'Potential V(x)', 'right_axis': True}
    wave_function_plot_abs_basic = {'color': 'purple', 'legend': 'Basic state wave function module'}
    wave_function_plot_abs_first = {'color': 'orange', 'legend': 'First state wave function module'}
    settings = {'name': 'zad_3_wave_100', 'title': 'Wave function of basic and first states for W = 2000 meV', 	'x_label': 'x [nm]', 'y_label': 'Wave functions module: $|\\Psi(x)|$', 'y_right_label': 'Potential V(x) [meV]', 'location': 'lower right'}
    Plot_1D([ wave_function_basic_x.real, wave_function_first_x.real, potential_x] ,
    [ abs(wave_function_basic_y), abs(wave_function_first_y), potential_y],
    [wave_function_plot_abs_basic, wave_function_plot_abs_first, potential_plot], settings, options['plot'], options['no_save']).plot()


    settings = {'name': 'zad_3_convergence_100', 'title': 'Documentation of the convergence of the solution for first state','x_label': 'Iteration', 'y_label': 'First state <E>', 'log': True, 'location': 'upper right'}
    Plot_1D([energy_x_basic, energy_x_first], [energy_y_basic, energy_y_first], [energy_plot_basic, energy_plot_first], settings, options['plot'], options['no_save']).plot()

    table = PrettyTable(['State', '<E> [meV]'])
    table.add_row(['Basic', energy_y_basic[-1]])
    table.add_row(['First', energy_y_first[-1]])
    print(table, '\n')

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--plot', action='store_true', help='Plot charts show')
    args.add_argument('--no_save', action='store_false', help='Plot charts save')
    args.add_argument('--debug', action='store_true', help='Debuging')
    options = vars(args.parse_args())

    parameters = numpy.arange(0.75, 1.2, 0.05)
    colors = ['red', 'coral', 'peru', 'olive', 'lime', 'forestgreen', 'blue', 'darkblue', 'black']
    colors = [{'basic_state': {'color': i}, 'first_state': {'color': i}} for i in colors]

    #zad_1(parameters, colors, options)
    #zad_2(options)
    zad_3(options)

