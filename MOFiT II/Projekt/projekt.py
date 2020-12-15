import sys
import time as TIME
import math
from scipy.linalg import eigh
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
        self._ax.get_figure().set_size_inches(10.5, 10.5, forward = True)

        if self._kwargs.setdefault('plot', True):
            plt.show()

    def _save_figure(self, name):

        self._ax.xaxis.label.set_size(18)
        self._ax.yaxis.label.set_size(18)
        self._ax.title.set_size(18)
        self._ax.tick_params(labelsize = 'large')
        self._ax.get_figure().set_size_inches(10.5, 10.5, forward = True)

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

class Covaled_Bond_Galerkin_Method:

    def __init__(self, **kwargs):
        self._alpha = kwargs['alpha']
        self._N = kwargs['N']
        self._a = kwargs['a']

    def _R_1(self, r):
        R = numpy.array([0,0,self._a])
        return numpy.linalg.norm(Atomic_Unit.nano2Bohr(r) - R)
    
    def _R_2(self, r):
        R = numpy.array([0,0,-self._a])
        return numpy.linalg.norm(Atomic_Unit.nano2Bohr(r) - R)

    def _ik(self, k):
        return 1 + (k-1)//self._N
    
    def _jk(self, k):
        return k - (self._ik(k) - 1) * self._N

    def _T_element(self, A, B, C, D):
        T = numpy.pi**(3/2)
        T1 = 3 * (A + B) * (C + D)
        T1 = T1 / ((A + B + C + D)**(5/2))
        T2 = 8 * self._a**2 * ((B*C - A*D)**2)
        T2 = T2 / ((A + B + C + D)**(7/2))
        T = T * (T1 - T2)
        T3 = -4 * self._a**2 * (A + C) * (B + D)
        T3 = T3 / (A + B + C + D) 
        T = T * numpy.exp(T3)
        return T

    def _V1_element(self, A, B, C, D):
        V = - numpy.pi**(3/2)
        V1 = 2 * self._a * (B + D) * numpy.sqrt(A + B + C + D)
        V = V / V1
        V2 = 2 * self._a * (B + D)
        V2 = V2 / numpy.sqrt(A + B + C + D)
        V = V * math.erf(V2)
        V3 = -4 * self._a**2 * (A + C) * (B + D)
        V3 = V3 / (A + B + C + D) 
        V = V * numpy.exp(V3)
        return V

    def _V2_element(self, A, B, C, D):
        V = - numpy.pi**(3/2)
        V1 = 2 * self._a * (A + C) * numpy.sqrt(A + B + C + D)
        V = V / V1
        V2 = 2 * self._a * (A + C)
        V2 = V2 / numpy.sqrt(A + B + C + D)
        V = V * math.erf(V2)
        V3 = -4 * self._a**2 * (A + C) * (B + D)
        V3 = V3 / (A + B + C + D) 
        V = V * numpy.exp(V3)
        return V

    def _H_matrix_element(self, k, l):
        A = self._alpha / (self._ik(k)**3)
        B = self._alpha / (self._jk(k)**3)
        C = self._alpha / (self._ik(l)**3)
        D = self._alpha / (self._jk(l)**3)
        T = self._T_element(A, B, C, D)
        V1 = self._V1_element(A, B, C, D)
        V2 = self._V2_element(A, B, C, D)
        return T + V1 + V2

    def _S_matrix_element(self, k, l):
        A = self._alpha / (self._ik(k)**3)
        B = self._alpha / (self._jk(k)**3)
        C = self._alpha / (self._ik(l)**3)
        D = self._alpha / (self._jk(l)**3)
        S = numpy.pi**(3/2)
        S1 = (A + B + C + D)**(3/2)
        S = S / S1
        S2 = -4 * self._a**2 * (A + C) * (B + D)
        S2 = S2 / (A + B + C + D)
        S = S * numpy.exp(S2)
        return S

    def H_matrix_fill(self):
        N2 = self._N**2
        H_kl = numpy.zeros((N2, N2))
        for k in range(N2):
            for l in range(N2):
                H_kl[k][l] = self._H_matrix_element(k + 1, l + 1)
        return H_kl

    def S_matrix_fill(self):
        N2 = self._N**2
        S_kl = numpy.zeros((N2, N2))
        for k in range(N2):
            for l in range(N2):
                S_kl[k][l] = self._S_matrix_element(k + 1, l + 1)
        return S_kl

    def Energy(self):
        eingenvalues, eingenvectors = eigh(self.H_matrix_fill(), self.S_matrix_fill())
        return (Atomic_Unit.Hartree2milieV(eingenvalues[0])/1000, numpy.rot90(eingenvectors)[::-1][0])

    def _fk(self, k, r1, r2):
        A = self._alpha / (self._ik(k)**3) * r1**2
        B = self._alpha / (self._jk(k)**3) * r2**2
        return numpy.exp(- A - B)

    def wave_function_value(self, r, vector):
        fi = 0
        N2 = self._N**2
        r1 = self._R_1(r)
        r2 = self._R_2(r)
        for k, ck in zip(range(N2), vector):
            fi += self._fk(k + 1, r1, r2) * 1
        return fi

def helper_zadanie_1(N, a, alpha_array):
    energy = [Covaled_Bond_Galerkin_Method(**{'alpha': alpha, 'N': N, 'a': a}).Energy()[0] for alpha in alpha_array]
    return {'x': alpha_array, 'y': energy, 'opt_energy': min(energy), 'opt_alpha': alpha_array[numpy.where(energy == min(energy))][0]}

def zadanie_1(options, name):
    print(f'\n{50*"-"} {name} {50*"-"}\n')
    start_time = TIME.time()
    resoults = []
    alpha_array = numpy.linspace(0.01, 15, 1000)
    for N in range(2,6):
        if options['debug']:
            print(f'Calculating for N = {N} in progress...')
        resoults.append({'N': N, **helper_zadanie_1(N, 1, alpha_array)})

    alpha_opt = [i['opt_alpha'] for i in resoults]
    energy_opt = [i['opt_energy'] for i in resoults]
    N = [i['N'] for i in resoults]
    X = [i['x'] for i in resoults]
    Y = [i['y'] for i in resoults]
    table = PrettyTable(['N', 'Optimal alpha', 'Energy [eV]'])
    for alpha, energy, n in zip(alpha_opt, energy_opt, N):
        table.add_row([n, alpha, '{:.4f}'.format(energy)])
    print('\n', table, '\n')
    print(f'{111*"_"}\n')

    colors = ['red', 'lime', 'blue', 'green', 'indigo']
    legends = []
    for i in N:
        legends.append({'color': colors[::-1][i-2], 'legend': 'N: {}'.format(i)})

    settings = {'name': 'wykres_alpha_najlepsze', 'title': 'Wykres zależności E w funkcji wartości parametru $\\alpha$ w zależności od N dla odległości a = 1 j. a.',
            'x_label': 'Wartość parametru $\\alpha$', 'y_label': 'Energia E [eV]', 'location': 'upper left'}
    Plot_1D(X, Y, [*legends], settings, options['plot'], options['no_save']).plot()

    return TIME.time() - start_time

def helper_zadanie_2(N, a_array, alpha):
    energy = [Atomic_Unit.Hartree2milieV(1/(2*a))/1000 + Covaled_Bond_Galerkin_Method(**{'alpha': alpha, 'N': N, 'a': a}).Energy()[0] for a in a_array]
    return {'x': a_array, 'y': energy, 'opt_energy': min(energy), 'opt_a': a_array[numpy.where(energy == min(energy))][0]}

def zadanie_2(options, name):
    print(f'\n{50*"-"} {name} {50*"-"}\n')
    start_time = TIME.time()
    resoults = []
    alpha_array = numpy.array([0.9103003003003003, 2.0206706706706705, 5.351781781781781, 10.108368368368367])
    a_array = numpy.linspace(2e-1, 4, 1000)
    for N in range(2, 6):
        if options['debug']:
            print(f'Calculating for N = {N} in progress...')
        resoults.append({'N': N, **helper_zadanie_2(N, a_array, alpha_array[N-2])})
    N = [i['N'] for i in resoults]
    X = [i['x'] for i in resoults]
    Y = [i['y'] for i in resoults]
   
    print(f'{111*"_"}\n')

    colors = ['red', 'lime', 'blue', 'green', 'indigo']
    legends = []
    for i in N:
        legends.append({'color': colors[::-1][i-2], 'legend': 'N: {}'.format(i)})

    settings = {'name': 'wykres_a_najlepsze', 'title': 'Wykres zależności $E_{tot}$ w funkcji odległości a [j. a.] od środka układu w zależności od N',
            'x_label': 'Wartość parametru a [j. a.]', 'y_label': 'Energia $E_{tot}$ [eV]', 'location': 'upper right'}
    Plot_1D(X, Y, [*legends], settings, options['plot'], options['no_save']).plot()

    return TIME.time() - start_time

def zadanie_3(options, name):
    print(f'\n{50*"-"} {name} {50*"-"}\n')
    start_time = TIME.time()
    resoults = []
    alpha = 5.351781781781781
    N = 4
    a_array = numpy.linspace(1e-1, 5, 1000)
    resoults = helper_zadanie_2(N, a_array, alpha)
    X = resoults['x']
    E_tot = resoults['y']
    E = E_tot - Atomic_Unit.Hartree2milieV(1/(2*a_array))/1000
    E_12 = Atomic_Unit.Hartree2milieV(1/(2*a_array))/1000
    E_teor = numpy.array([Atomic_Unit.Hartree2milieV(-0.6026)/1000 for _ in X])
    label = '$E_{teor}$ = ' + '{:.4f} [eV]'.format(Atomic_Unit.Hartree2milieV(-0.6026)/1000)
    legends = [
        {'color': 'red', 'legend': '$E_{tot}$'},
        {'color': 'blue', 'legend': '$E$'},
        {'color': 'green', 'legend': '$\\frac{1}{2a}$  [eV]'},
        {'color': 'black', 'legend': label}
    ]
    print('Wartość minimalna E_tot {:.4f} [eV] dla a = {:.3f} [j. a.]\n'.format(min(E_tot), X[numpy.where(E_tot == min(E_tot))][0]))
    settings = {'name': 'wykres_E_tot', 'title': 'Wykres zależności $E$ w funkcji odległości a [j. a.] od środka układu dla N = 4',
            'x_label': 'Wartość parametru a [j. a.]', 'y_label': 'Energia $E$ [eV]', 'location': 'upper right'}
    Plot_1D([X, X, X, X], [E_tot, E, E_12, E_teor], [*legends], settings, options['plot'], options['no_save']).plot()
    return TIME.time() - start_time
    
def zadanie_4(options, name):
    print(f'\n{50*"-"} {name} {50*"-"}\n')
    start_time = TIME.time()
    resoults = []
    alpha = 5.351781781781781
    N = 4
    a_array = numpy.linspace(1e-1, 2.5, 1000)
    resoults = helper_zadanie_2(N, 2*a_array, alpha)
    X = resoults['x']
    E_tot = resoults['y']
    E_bind = - numpy.array(E_tot) - Atomic_Unit.Hartree2milieV(1/2)/1000 
    E_teor = numpy.array([2.793 for _ in X])

    l = [(a, Eb) for a, Eb, Et in zip(X, E_bind, E_tot) if Eb >= Et ]
    print('Wartość minimalna 2a {:.3f} (E_bind = {:.4f} [ev]). Wartość maksymalna 2a {:.3f} (E_bind = {:.4f} [ev])\n'.format(l[0][0], l[0][1], l[-1][0], l[-1][1]))

    legends = [
        {'color': 'red', 'legend': '$E_{tot}$'},
        {'color': 'blue', 'legend': '$E_{bind}$'},
        {'color': 'black', 'legend': '$E_{bind}^{teor}$ = 2.793 [eV]'}
    ]
    settings = {'name': 'wykres_E_bind', 'title': 'Wykres zależności $E$ w funkcji odległości 2a [j. a.] dla N = 4',
            'x_label': 'Wartość parametru 2a', 'y_label': 'Energia $E$ [eV]', 'location': 'upper right'}
    Plot_1D([X, X, X], [E_tot, E_bind, E_teor], [*legends], settings, options['plot'], options['no_save']).plot()
    return TIME.time() - start_time

def zadanie_5(options, name):
    print(f'\n{50*"-"} {name} {50*"-"}\n')
    start_time = TIME.time()
    Object = Covaled_Bond_Galerkin_Method(**{'N': 4, 'a': 1, 'alpha': 5.351781781781781})
    E = Object.Energy()
    X = numpy.arange(-0.2, 0.2, 0.001)
    Z = numpy.arange(-0.2, 0.2, 0.001)
    fi = []
    for i, x in enumerate(X):
        fi.append([])
        for j, z in enumerate(Z):
            r = numpy.array([x, 0, z])
            q = (Object.wave_function_value(r, E[1]))
            fi[i].append(q)

    options = {'name' : 'rozklad_ladunku', 'title' : 'Rozkład ładunku wyrażony jako funkcja $\\varphi(\\vec{r} = (x,0,z))$','xlabel': 'Oś x [nm]', 'ylabel': 'Oś z [nm]' ,
            'no_save': True, 'plot': False,
            'colorbar_title': '$\\varphi(\\vec{r})$', 'extent': [-0.2, 0.2, -0.2, 0.2], 'heatmap_color': 'BuPu'}

    Plot_Imshow(fi, **options).plot()



    return TIME.time() - start_time

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--plot', action='store_true', help='Plot charts show')
    args.add_argument('--no_save', action='store_false', help='Plot charts save')
    args.add_argument('--debug', action='store_true', help='Debuging')
    options = vars(args.parse_args())

    time_table = []
    time_table.append(('1', zadanie_1(options, 'ZADANIE 1')))
    time_table.append(('2', zadanie_2(options, 'ZADANIE 2')))
    time_table.append(('3', zadanie_3(options, 'ZADANIE 3')))
    time_table.append(('4', zadanie_4(options, 'ZADANIE 4')))
    time_table.append(('5', zadanie_5(options, 'ZADANIE 5')))

    table = PrettyTable(['Task', 'Time [s]'])
    for i, value in time_table:
        table.add_row([i, value])
    print(table, '\n')
