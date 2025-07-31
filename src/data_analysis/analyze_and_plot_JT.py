import os
import sys
import git
import math
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from scipy import optimize, interpolate
from numpy.polynomial.chebyshev import Chebyshev
from FluidProps import FluidProps
from get_git_root import get_git_root

os_idx = ["linux", "win32", "darwin"].index(sys.platform)
root = ["/etc/REFPROP", "c:/Program Files (x86)/REFPROP", "/Users/jtkaczuk/codes/REFPROP"][os_idx]
if not os.path.exists(root):
    print(f"{root} does not exist. Set existing path to REFPROP.")
    raise

path_panorama = os.path.join(path, "PanoramaDataAcquisition")
path_exp = os.path.join(path_panorama, "JT_measurements")
path_Ar = os.path.join(path, "Ar_CSV")
path_He = os.path.join(path, "He_CSV")
path_N2 = os.path.join(path, "N2_CSV")
path_HeNe = os.path.join(path, "HeNe_CSV")
path_HeN2 = os.path.join(path, "HeN2_CSV")


def set_mpl():
    fontsize = 16
    mpl.rcParams.update({"font.family": "Times New Roman",
                         "mathtext.fontset": "dejavuserif",
                         "font.size": fontsize,
                         "axes.labelsize": fontsize,
                         "axes.titlesize" : fontsize,
                         "legend.fontsize": fontsize,
                         "xtick.top": True,
                         "xtick.bottom": True,
                         "ytick.left": True,
                         "ytick.right": True,
                         "xtick.direction": "in",
                         "ytick.direction": "in",
                         "xtick.major.pad": 5,
                         "ytick.major.pad": 5})


class CernoxCal(object):
    serNum = str

    def __init__(self, serNum):
        assert serNum in ['X93303','X115143','X115888']
        self.serNum = serNum
        self.coefficients = self.set_coefficients()
        self.tbl = self.set_tbl()
        self.dat = self.set_dat()

    def set_coefficients(self):
        fpath = os.path.join(self.get_git_root(),
                             "data",
                             "cernox_calibration_data",
                             self.serNum,
                             f"{self.serNum}.cof")
        d = {}
        with open(fpath, 'r') as fh:
            coefs = []
            limits = ''
            while True:
                # read the calibration file to get calibration curve coefficient
                line = fh.readline()
                if not line:
                    d[limits] = {'Zlower': Zlower, 'Zupper': Zupper, 'coefs': coefs}
                    break
                l = line.strip().split()
                if l[0] == 'Lower' and l[1] == 'Resist.':
                    limits = l[-1]
                elif l[0] == 'Upper' and l[1] == 'Resist.':
                    limits = limits + ',' + l[-1]
                elif l[0] == 'Zlower':
                    Zlower = float(l[-1])
                elif l[0] == 'Zupper':
                    Zupper = float(l[-1])
                elif list(l[0])[0] == 'C' and list(l[0])[1] == '(':
                    coefs.append(float(l[-1]))
                elif coefs and list(l[0])[0] != 'C':
                    d[limits] = {'Zlower': Zlower,
                                 'Zupper': Zupper,
                                 'coefs': coefs}
                    coefs = []
        return d

    def set_tbl(self):
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def set_dat(self):
        fpath = os.path.join(self.get_git_root(),
                     "data",
                     "cernox-calibration-data",
                     self.serNum,
                     f"{self.serNum}.tbl")
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def get_T_from_coefs(self, R):
        """
        calculate temperature in K from:
          - resistance in Ohm
          - coefficients from .cof file
        """
        T = 0
        for key, value in self.coefficients.items():
            L, U = key.split(',')
            L, U = float(L), float(U)
            if R >= L and R < U:
                k = ((np.log10(R) - value['Zlower']) -
                    (value['Zupper'] - np.log10(R))) / (
                    value['Zupper'] - value['Zlower'])
                for i, c in enumerate(value['coefs']):
                    T += c * np.cos(i * np.arccos(k))
        return T

    def get_T_from_tbl(self, R):
        """
        calculate temperature in K from resistance in Ohm
        by interpolating values from tab file
        """
        f = interpolate.interp1d(self.tbl[('Resistance', '(Ohms)')],
                                 self.tbl[('Temp.', '(Kelvin)')])
        return f(R)

    def get_interp_abs_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Tt = self.get_T_from_tbl(R)
        return Tc - Tt

    def get_interp_rel_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Tt = self.get_T_from_tbl(R)
        return (Tc - Tt) / Tt

    def get_deltaT_dat_cof(self):
        T, DT = [], []
        for row in self.dat.iterrows():
            Tdat = row[1][('Temperature','(Kelvin)')]
            T.append(Tdat)
            DT.append(Tdat - self.get_T_from_coefs(row[1][('Resistance','(Ohms)')]))
        return T, DT

    def get_polynomial_fit_uncert(self, T):
        """ Chebyshev polynomial fitting standard expanded uncertainty in K """
        d = {
            'X93303': [{'Tmin': 1.4 , 'Tmax': 14.1 , 'N': 31, 'n': 9, 'DTrms': 0.93},
                      {'Tmin': 14.1, 'Tmax': 80.0 , 'N': 31, 'n': 7, 'DTrms': 1.60},
                      {'Tmin': 80.0, 'Tmax': 325.0, 'N': 32, 'n': 8, 'DTrms': 5.84}],
            'X115143': [{'Tmin': 20.0, 'Tmax': 95.0, 'N': 28, 'n': 7, 'DTrms': 0.91},
                       {'Tmin': 95.0, 'Tmax': 325.0, 'N': 29, 'n': 9, 'DTrms': 3.45}],
            'X115888': [{'Tmin': 20.0, 'Tmax': 95.3, 'N': 28, 'n': 7, 'DTrms': 1.17},
                       {'Tmin': 95.3, 'Tmax': 325.0, 'N': 29, 'n': 8, 'DTrms': 4.06}]
        }
        for c in d[self.serNum]:
            if T >= c['Tmin'] and T < c['Tmax']:
                std = (c['N'] / (c['N'] - c['n']) * c['DTrms']**2)**0.5
                N = c['N']
                break
            else:
                std = 0
                N = 1
        return std / N * 1e3


class StdTempUncertainty(object):
    sensor = str

    def __init__(self, sensor):
        self.cernox = CernoxCal(sensor)

    def get_R_from_T(self, T):
        dfx0 = self.cernox.dat.iloc[(self.cernox.dat[('Temperature', '(Kelvin)')] - T).abs().argsort()[:1]]
        res = optimize.minimize(fun = lambda x: np.abs(self.cernox.get_T_from_coefs(x) - T),
                                x0 = dfx0[('Resistance', '(Ohms)')])
        return res.x[0]

    def get_cernox_uncertainty(self, T):
        """ temperature measurement uncertainty of the cernox sensor in K (between 20 K and 300 K) """
        assert T >= 20 and T <= 300
        return (-5e-05 * T**2 + 0.1281 * T + 6.4603) / 1e3

    def get_cabtr_uncertainty(self, T):
        """ temperature uncertainty of data acquisition module in K (between 20 K and 300 K) """
        assert T >= 20 and T <= 300
        return (2e-05 * T**2 + 0.0695 * T - 0.1001) / 1e3

    def __call__(self, T):
        cabtr = self.get_cabtr_uncertainty(T)
        cernox = self.get_cernox_uncertainty(T)
        poly = self.cernox.get_polynomial_fit_uncert(T)
        return (cabtr**2 + cernox**2 + poly**2)**0.5


class Isenthalpic(object):
    fluid = str or tuple
    Dint = float
    k = int
    pacc = float
    pmax = float
    
    def __init__(self,
                 fluid,
                 Dint = 0.365e-3,
                 k = 1.96,
                 pacc = 1e-4,
                 pmax = 13.7):
        self.fluid = FluidProps(fluid)
        self.pacc = pacc
        self.pmax = pmax
        self.Dint = Dint
        self.Aflow = np.pi * Dint**2 / 4
        self.k = k

    def get_u(self, m, p, T):
        """ return fluid velocity in m/s """
        return m / self.fluid.get_rhomass(p, T) / self.Aflow
        
    def get_DT2h(self, m, p1, p2, T1, T2h):
        """ return temeprature difference between isenthalpic and non-isenthalpic conditions """
        T2 = T2h - (self.get_u(m=m, p=p2, T=T2h)**2
            - self.get_u(m=m, p=p1, T=T1)**2) / (2 * self.fluid.get_cp(p=p2, T=T2h))
        return T2h - T2

    def get_T2_corrected(self, p1, p2, T1, T2, x=1.0):
        fluidList = ['helium','nitrogen','argon','helium;neon','helium;nitrogen']
        idx = fluidList.index(self.fluid.lower())
        path = [path_He,path_N2,path_Ar,path_HeNe,path_HeN2][idx]
        fldName = ['He','N2','Ar','HeNe','HeN2'][idx]
        if idx in [0, 1, 2]:
            fileName = '{}_D0.365_L6_p{:.1f}_T{:.1f}.csv'.format(fldName, p1, T1)
        else:
            fileName = '{}_D0.365_L6_p{:.1f}_T{:.1f}_{}He.csv'.format(fldName, p1, T1, round(x, 1))
        df = pd.read_csv(os.path.join(path, fileName), sep=',', engine='python')
        interp = interpolate.interp1d(df['p2/MPa'], df['m/kg_per_s'])
        if p2 < 0.5:
            return T2
        m = float(interp(p2))
        u1 = self.get_u(m=m, p=p1, T=T1)
        u2 = self.get_u(m=m, p=p2, T=T2)
        T2_corrected = T2 - (u2**2 - u1**2) / (2 * self.get_cp(p=p2, T=T2))
        if T2_corrected <= self.get_saturation_T(p2):
            return T2
        else:
            return T2_corrected
       
    def get_JT_err(self, p1, p2, T1):
        TUncertIn = TempUncertainty('X115143')
        TUncertOut = TempUncertainty('X93303')

        T2h = self.getT(p=p2, h=self.geth(p=p1, T=T1))
        DTnonh = self.get_DT2h(m=m, p1=p1, p2=p2, T1=T1, T2h=T2h)
        T2 = T2h - DTnonh

        U_P = self.pacc * self.pmax # MPa

        U_Tin = TUncertIn(T1) # K
        U_Tout = TUncertOut(T2) # K

        return ((U_Tin**2 + (U_Tout + self.k * DTnonh)**2) / (T1 - T2)**2
            + 2 * U_P**2 / (p1 - p2)**2)**0.5 * 100

    def get_Ur_muJT(self, pin, pout, Tin, Tout):
        TUncertIn = TempUncertainty('X115143')
        TUncertOut = TempUncertainty('X93303')

        U_P = self.pacc * self.pmax # MPa

        U_Tin = TUncertIn(Tin) # K
        U_Tout = TUncertOut(Tout) # K

        return ((U_Tin**2 + U_Tout**2) / (Tin - Tout)**2
            + 2 * U_P**2 / (pin - pout)**2)**0.5 * 100


def plot_isenthalps_and_JT_coefs(fluid,
                                 fnames,
                                 order=None,
                                 errorBars=False,
                                 uncertaintyAnalysis='conventional',
                                 saturationLine=False,
                                 **kwargs):
    fig, [ax0, ax1, ax2] = plt.subplots(3, 1,
                                        figsize = (3.2, 7.5),
                                        sharex = True,
                                        gridspec_kw = {'height_ratios':[1.8, 1.2, 0.8]})
    ax2.hlines(0, 0, 15, linewidth=1., color='k', zorder=-100)
    colors = ('blue', 'red', 'green', 'purple', 'orange',
              'deeppink', 'yellowgreen', 'navy', 'coral', 'tomato')
    isenthalp = Isenthalpic(fluid)
    
    res_d = {}
    uJT_list, pin_list, Tin_list, p_list, T_list = [], [], [], [], []
    x1_list, x2_list = [], []

    if order is None:
        orders = [None]*len(fnames)
    elif (type(order) is int) or (type(order) is float):
        orders = [order]*len(fnames)
    else:
        orders = order

    for i, fname, color, order in zip(range(len(fnames)), fnames, colors, orders):
        df = pd.read_csv(fname, sep=',')
        if type(fluid) is tuple:
            xHe = np.mean(df['x(He)'])
            isenthalp.set_composition_from_1st_fraction(xHe)
        else:
            xHe = 1.0

        if order is None:
            deg_d = {}
            for deg in range(len(df['PT102/MPa'])):
                newSeries = Chebyshev.fit(x = df['PT102/MPa'],
                                           y = df['TT102/K'],
                                           deg = deg)
                fit_ = Chebyshev(newSeries.convert().coef)
                newSeriesDeriv = fit_.deriv(m=1)
                derivFit_ = Chebyshev(newSeriesDeriv.convert().coef)

                err_JT = []
                for _, row in df.iterrows():
                    uJT_meas = derivFit_(row['PT102/MPa'])
                    uJT_calc = isenthalp.get_JT_coefficient(p = row['PT102/MPa'],
                                                           T = row['TT102/K'])
                    err_JT.append(abs((uJT_calc - uJT_meas) / uJT_meas))
                deg_d[deg] = sum(err_JT)
            deg = int(min(deg_d, key=deg_d.get))
        else:
            deg = order
        print('Polynomial degree: {}, number of points: {}'.format(deg, len(df['PT102/MPa'])))

        # calculations necessary for plotting
        _h = np.array([])
        for _p, _T in zip(df['PT102/MPa'], df['TT102/K']):
            _h = np.append(_h, isenthalp.get_hmass(p=_p, T=_T))
        df['h2/(J/kg)'] = _h
        p = np.arange(df['PT102/MPa'].min(), df['PT102/MPa'].max()+0.05, 0.05)
        hmins = []
        for _h in df['h2/(J/kg)']:
            T = np.array([isenthalp.get_T(p=_, h=_h) for _ in p])
            hmins.append(sum([min(abs(T - _T)) for _T in df['TT102/K']]))
        ind = hmins.index(min(hmins))
        T = np.array([isenthalp.get_T(p=_, h=df['h2/(J/kg)'][ind]) for _ in p])

        newSeries = Chebyshev.fit(x = df['PT102/MPa'],
                                   y = df['TT102/K'],
                                   deg = deg)
        fit_ = Chebyshev(newSeries.convert().coef)
        newSeriesDeriv = fit_.deriv(m=1)
        derivFit_ = Chebyshev(newSeriesDeriv.convert().coef)

        if type(fluid) is str:
            print('pin/MPa, Tin/K, p/MPa, T/K, uJT_meas/(K/MPa), uJT_calc/(K/MPa), error')
        else:
            print('pin/MPa, Tin/K, p/MPa, T/K, x1, x2, uJT_meas/(K/MPa), uJT_calc/(K/MPa), error')
        DT, err_JT, _uJT_meas = [], [], []
        
        # calculate errors
        for j, row in enumerate(df.iterrows()):
            T_meas = row[1]['TT102/K']
            T_calc = isenthalp.get_T(p = row[1]['PT102/MPa'],
                                     h = row[1]['TT102/K'])
            DT.append(T_meas - T_calc)
            uJT_meas = derivFit_(row[1]['PT102/MPa'])
            _uJT_meas.append(uJT_meas)
            uJT_calc = isenthalp.get_JT_coefficient(p = row[1]['PT102/MPa'],
                                                   T = row[1]['TT102/K'])
            err_JT.append((uJT_calc - uJT_meas) / uJT_meas * 100)
            if j != 0 and j != (df.shape[0]-1):
                pin_list.append(row[1]['PT101/MPa'])
                Tin_list.append(row[1]['TT101/K'])
                p_list.append(row[1]['PT102/MPa'])
                T_list.append(row[1]['TT102/K'])
                x1_list.append(row[1]['x(He)'])
                x2_list.append(1-row[1]['x(He)'])
                uJT_list.append(uJT_meas)
                if type(fluid) is str:
                    print('{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}'.format(
                          df['PT101/MPa'][0], df['TT101/K'][0],
                          row[1]['PT102/MPa'], row[1]['TT102/K'],
                          uJT_meas, uJT_calc, err_JT[-1]/100))
                else:
                    print('{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.4f}, {:.4f}, {:.6f}, {:.6f}, {:.6f}'.format(
                          df['PT101/MPa'][0], df['TT101/K'][0],
                          row[1]['PT102/MPa'], row[1]['TT102/K'],
                          df['x(He)'][0], 1-df['x(He)'][0],
                          uJT_meas, uJT_calc, err_JT[-1]/100))
        print('\n')

        # plot isenthalpic line, measurement points, and fitted line
        if i == 0:
            eos_label, m_label, c_label = 'EOS', 'measurements', 'fit'
        else:
            eos_label, m_label, c_label = None, None, None
        ax0.plot(df['PT102/MPa'][1:-1], df['TT102/K'][1:-1],
                 linestyle = '', color = color,
                 marker = 'o', ms = 8, mew = 1.1,
                 mfc = 'w', label = m_label)
        ax0.plot(df['PT102/MPa'][0], df['TT102/K'][0], 
                 df['PT102/MPa'].values[-1], df['TT102/K'].values[-1],
                 linestyle = '', color = color,
                 marker = 'o', ms = 8, mew = 1.1,
                 mfc = color, label = None)
        ax0.plot(p, T, linewidth=1., color = color, label = eos_label)
        ax0.plot(p, fit_(p),
                 linewidth = 1., linestyle = ':',
                 color = color, label = c_label)
        _p = np.arange(0.1, 15.1, 0.1)
        _T = np.array([isenthalp.get_T(p=_, h=df['h2/(J/kg)'][ind]) for _ in _p])
        _u = np.array([isenthalp.get_JT_coefficient(p=__p, T=__T) for __p, __T in zip(_p, _T)])
        
        # remove plotting errors when getting too close to the the saturation line
        if max(_u[1:] - _u[:-1]) < 0.1:
            ax1.plot(_p, _u, linewidth = 1., color = color)
        else:
            _p = np.arange(df['PT102/MPa'].min(), 15.05, 0.05)
            _T = np.array([isenthalp.get_T(p=_, h=df['h2/(J/kg)'][ind]) for _ in _p])
            _u = np.array([isenthalp.get_JT_coefficient(p=__p, T=__T) for __p, __T in zip(_p, _T)])
            ax1.plot(_p, _u, linewidth = 1., color = color, zorder=10)
        
        # plot dotted line in the Joule-Thomson coefficient figure
        intMin = min(df['PT102/MPa'][1:-1])
        intMax = max(df['PT102/MPa'][1:-1])
        _p = np.arange(intMin, intMax+0.1, 0.1)
        ax1.plot(_p, derivFit_(_p), linewidth = 1.,
                 linestyle = ':', color = color, zorder=10)
        ax1.plot(df['PT102/MPa'][1:-1], _uJT_meas[1:-1],
                 linewidth = 1., linestyle = '',
                 color = color, marker = 'o', ms = 8,
                 mfc = 'w', mew =1.1, label = None, zorder=-10)
        
        # plot error bars
        if errorBars:
            if type(fluid) is str:
                df_res = pd.read_csv(os.path.join(path_panorama, 'MEASUREMENT_pures.csv'))
                df_res = df_res[df_res['fluid'] == fluid]
            else:
                df_res = pd.read_csv(os.path.join(path_panorama, 'MEASUREMENT_mixtures.csv'))
                df_res = df_res[df_res['fluid'] == '{}-{}'.format(*fluid)]
            df_plot = df_res[round(df_res['pin/MPa']) == round(df['PT101/MPa'][0])]
            df_plot = df_plot[round(df_plot['Tin/K']) == round(df['TT101/K'][0])]
            if type(fluid) is tuple:
                xRounded = float(fname.rstrip('.csv').split('_')[-1].replace('x',''))
                df_plot = df_plot[df_plot['x1_csvName'] == xRounded]
            
            if uncertaintyAnalysis == 'conventional':
                ax2.errorbar(df['PT102/MPa'][1:-1], err_JT[1:-1],
                    yerr = df_plot['combRelU'] * 100,
                    linestyle = '', color = color, mfc = 'w',
                    ms = 8, mew=.9, marker = 'o', ecolor = 'k',
                    elinewidth = 1.1, capsize = 2, barsabove = True)
            elif uncertaintyAnalysis == 'MonteCarlo':
                ax2.errorbar(df['PT102/MPa'][1:-1], err_JT[1:-1],
                    yerr = df_plot['MonteCarloUncert'] * 100,
                    linestyle = '', color = color, mfc = 'w',
                    ms = 8, mew=.9, marker = 'o', ecolor = 'k',
                    elinewidth = 1.1, capsize = 2, barsabove = True)
            elif uncertaintyAnalysis == 'both':
                yerr_conv =  df_plot['combRelU'] * 100
                yerr_MC = df_plot['MonteCarloUncert'] * 100
                ax2.errorbar(df['PT102/MPa'][1:-1], err_JT[1:-1],
                    yerr = yerr_MC, linestyle = '', color = color,
                    mfc = 'w', ms = 8, mew=.9, marker = 'o', ecolor = 'gray',
                    elinewidth = 1.1, capsize = 2, barsabove = True)
                ax2.errorbar(df['PT102/MPa'][1:-1], err_JT[1:-1],
                    yerr = yerr_conv, linestyle = '', marker = '', ecolor = 'k',
                    elinewidth = 1.1, capsize = 2, barsabove = True)
        
        else:
            ax2.plot(df['PT102/MPa'][1:-1], err_JT[1:-1],
                 linestyle = '', marker = 'o',
                 ms = 5.5, color = color, mfc = 'w')
        
        # plot saturation line
        if saturationLine:
            pc, Tc = isenthalp.get_critical_p(), isenthalp.get_critical_T()
            psat = np.arange(0.1, round(pc)+1, 0.005)
            Tsat = np.array([isenthalp.get_saturation_T(_) for _ in psat])
            ax0.plot(psat[Tsat > 0], Tsat[Tsat > 0], color='k', linewidth=1., zorder=-20)
            ax0.plot(pc, Tc, color='k', marker='o', ms=4.5)

    # set limits on axes
    if 'xlim' in kwargs:
        ax0.set_xlim(kwargs['xlim'])
    if 'ylim0' in kwargs:
        ax0.set_ylim(kwargs['ylim0'])
    if 'xticks' in kwargs:
        ax0.set_xticks(kwargs['xticks'])
    if 'yticks0' in kwargs:
        ax0.set_yticks(kwargs['yticks0'])
    ax0.set_ylabel(r'$T~/~{\rm K}$', labelpad=6) # 10
    if 'legend' in kwargs:
        if kwargs['legend'] is True:
            ax0.legend(frameon=False,
                       handletextpad=0.4,
                       bbox_to_anchor=(-0.02,-0.06,1,1))
                       # loc='lower center')
    if type(fluid) is tuple:
        ax0.set_title(r'$x_{{\rm He}} = {:.4f}$'.format(xHe))
    ax1.set_ylabel(r'$\mu_{\rm JT}~/~\left({\rm K~MPa^{-1}}\right)$', labelpad=10) # 0
    ax2.set_xlabel(r'$p~/~{\rm MPa}$')
    ax2.set_ylabel(r'$100~\Delta \mu_{\rm JT}~/~\mu_{\rm JT}$', labelpad=-4) # 5
    if 'ylim1' in kwargs:
        ax1.set_ylim(kwargs['ylim1'])
    if 'yticks1' in kwargs:
        ax1.set_yticks(kwargs['yticks1'])
    if 'ylim2' in kwargs:
        ax2.set_ylim(kwargs['ylim2'])
    if 'yticks2' in kwargs:
        ax2.set_yticks(kwargs['yticks2'])
    fig.tight_layout(pad=0.2)

    res_d['p_in/MPa'] = pin_list
    res_d['T_in/MPa'] = Tin_list
    res_d['p/MPa'] = p_list
    res_d['T/K'] = T_list
    if type(fluid) is tuple:
        res_d['x1'] = x1_list
        res_d['x2'] = x2_list
    res_d['uJT/(K/MPa)'] = uJT_list

    return res_d


def validate_nitrogen(uncertaintyAnalysis='conventional'):
    """ Method validation with nitrogen """
    print('Nitrogen')
    # fname1 = os.path.join(path_exp,'20201130_N2_140K_5MPa.csv')
    fname2 = os.path.join(path_exp,'20210304_N2_hysteresis_01.csv')
    fname3 = os.path.join(path_exp,'20201214_N2_160K_5MPa.csv')
    fname4 = os.path.join(path_exp,'20201217_N2_160K_6MPa.csv')
    fname5 = os.path.join(path_exp,'20201217_N2_160K_7MPa.csv')
    # fname5 = os.path.join(path_exp,'20210125_N2_160K_7MPa.csv')
    fname6 = os.path.join(path_exp,'20201127_N2_160K_9MPa.csv')
    fname7 = os.path.join(path_exp,'20201216_N2_160K_12MPa.csv')
    kwargs = {'xlim': (0, 12),
          'ylim0': (80, 165),
          'ylim1': (0, 20),
          'ylim2': (-3, 3),
          'xticks': np.arange(0, 13, 3),
          'yticks0': np.arange(80, 165, 20),
          'yticks1': np.arange(0, 25, 5),
          'yticks2': np.arange(-3, 4, 1.5),
          'legend': True}
    res = plot_isenthalps_and_JT_coefs(
        fluid = 'Nitrogen',
        fnames = [fname2, fname3, fname4, fname5, fname6, fname7],
        order = [None, None, 4, None, None, None],
        errorBars = True,
        saturationLine = False,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs
        )
    df = pd.DataFrame.from_dict(res)
    df.to_csv('JT_res_N2.csv', sep=r',', index=False)


def validate_argon(uncertaintyAnalysis='conventional'):
    """ Method validation with argon """
    print('Argon')
    fname1 = os.path.join(path_exp,'20201215_Ar_180K_12MPa.csv')
    fname2 = os.path.join(path_exp,'20201215_Ar_180K_5MPa.csv')
    kwargs = {'xlim': (0, 15),
              'ylim0': (140, 190),
              'ylim1': (0, 15),
              'ylim2': (-3, 3),
              'xticks': np.arange(0, 16, 5),
              'yticks0': np.arange(140, 200, 10),
              'yticks1': np.arange(0, 16, 5),
              'yticks2': np.arange(-3, 4, 1.5)}
    res = plot_isenthalps_and_JT_coefs(
        fluid = 'Argon',
        fnames = [fname1, fname2],
        errorBars = True,
        saturationLine = False,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs
        )
    df = pd.DataFrame.from_dict(res)
    df.to_csv('JT_res_Ar.csv', sep=r',', index=False)

def validate_helium(uncertaintyAnalysis='conventional'):
    """ method validation with helium """
    print('Helium')
    fname1 = os.path.join(path_exp,'20210112_He_65K_7MPa.csv')
    fname2 = os.path.join(path_exp,'20210113_He_65K_5MPa.csv')
    fname3 = os.path.join(path_exp,'20210126_He_140K_10MPa.csv')
    kwargs = {'xlim': (0, 12),
              'ylim0': (60, 160),
              'ylim1': (-0.8, 0),
              'ylim2': (-20, 20),
              'xticks': np.arange(0, 13, 3),
              'yticks0': np.arange(60, 165, 20),
              'yticks1': np.arange(-0.8, 0.1, 0.2),
              'yticks2': np.arange(-20, 30, 10)}
    res = plot_isenthalps_and_JT_coefs(
        fluid = 'Helium',
        fnames = [fname1, fname2, fname3],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs
        )
    df = pd.DataFrame.from_dict(res)
    df.to_csv('JT_res_He.csv', sep=r',', index=False)


def JT_helium_nitrogen(uncertaintyAnalysis='conventional'):
    """ plot isenthalps and measurements deviation from equation """
    print('Helium-Nitrogen')
    fname1 = os.path.join(path_exp,'20201209_HeN2_160K_8MPa_x0.145.csv')
    fname2 = os.path.join(path_exp,'20201210_HeN2_160K_5MPa_x0.145.csv')
    fname3 = os.path.join(path_exp,'20201210_HeN2_140K_5MPa_x0.145.csv')
    fname4 = os.path.join(path_exp,'20201211_HeN2_140K_5MPa_x0.503.csv')
    fname5 = os.path.join(path_exp,'20201211_HeN2_140K_8MPa_x0.503.csv')
    kwargs1 = {'xlim': (0, 10),
               'ylim0': (110, 165),
               'ylim1': (0, 10),
               'ylim2': (-4, 4),
               'xticks': np.arange(0, 11, 2),
               'yticks2': np.arange(-4, 5, 2)}
    kwargs2 = {'xlim': (0, 10),
               'ylim0': (110, 165),
               'ylim1': (0, 5),
               'ylim2': (-20, 20),
               'xticks': np.arange(0, 11, 2),
               'yticks2': np.arange(-20, 30, 10)}
    res1 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Nitrogen'),
        fnames = [fname1, fname2, fname3],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs1
        )
    res2 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Nitrogen'),
        fnames = [fname4, fname5],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs2
        )
    df = pd.DataFrame.from_dict(res1)
    df.to_csv('JT_res_HeN2_01.csv', sep=r',', index=False)
    df = pd.DataFrame.from_dict(res2)
    df.to_csv('JT_res_HeN2_02.csv', sep=r',', index=False)


def JT_helium_neon(uncertaintyAnalysis='conventional'):
    print('Helium-Neon')
    fname1 = os.path.join(path_exp,'20210114_HeNe_65K_10MPa_x0.46.csv')
    fname2 = os.path.join(path_exp,'20210114_HeNe_80K_10MPa_x0.46.csv')
    fname3 = os.path.join(path_exp,'20210115_HeNe_65K_8MPa_x0.46.csv')
    fname4 = os.path.join(path_exp,'20210115_HeNe_80K_5MPa_x0.46.csv')
    fname5 = os.path.join(path_exp,'20210115_HeNe_80K_7MPa_x0.46.csv')
    fname6 = os.path.join(path_exp,'20210127_HeNe_65K_10MPa_x0.32.csv')
    fname7 = os.path.join(path_exp,'20210127_HeNe_80K_8MPa_x0.32.csv')
    fname8 = os.path.join(path_exp,'20210128_HeNe_65K_5MPa_x0.32.csv')
    fname9 = os.path.join(path_exp,'20210128_HeNe_65K_5MPa_x0.4.csv')
    fname10 = os.path.join(path_exp,'20210128_HeNe_65K_7MPa_x0.4.csv')
    fname11 = os.path.join(path_exp,'20210303_HeNe_65K_5MPa_x0.2.csv')
    fname12 = os.path.join(path_exp, '20210303_HeNe_65K_7MPa_x0.2.csv')
    kwargs = {'xlim': (0, 10),
              'ylim0': (50, 85),
              'ylim1': (0, 4),
              'ylim2': (-10, 10),
              'xticks': np.arange(0, 11, 2),
              'yticks1': np.arange(0, 5, 1),
              'yticks2': np.arange(-10, 11, 5)}
    kwargs_leg = kwargs.copy()
    kwargs_leg['legend'] = True
    res1 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Neon'),
        fnames = [fname1, fname2, fname3, fname4, fname5],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs)
    res2 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Neon'),
        fnames = [fname6, fname7, fname8],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs)
    res3 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Neon'),
        fnames = [fname9, fname10],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs)
    res4 = plot_isenthalps_and_JT_coefs(
        fluid = ('Helium', 'Neon'),
        fnames = [fname11, fname12],
        order = None,
        errorBars = True,
        uncertaintyAnalysis = uncertaintyAnalysis,
        **kwargs_leg)
    df = pd.DataFrame.from_dict(res1)
    df.to_csv('JT_res_HeNe_01.csv', sep=r',', index=False)
    df = pd.DataFrame.from_dict(res2)
    df.to_csv('JT_res_HeNe_02.csv', sep=r',', index=False)
    df = pd.DataFrame.from_dict(res3)
    df.to_csv('JT_res_HeNe_03.csv', sep=r',', index=False)
    df = pd.DataFrame.from_dict(res4)
    df.to_csv('JT_res_HeNe_04.csv', sep=r',', index=False)

if __name__ == '__main__':
    set_mpl()
    # uncertaintyAnalysis = 'MonteCarlo'
    # uncertaintyAnalysis = 'conventional'
    uncertaintyAnalysis = 'both'

    # validate_nitrogen(uncertaintyAnalysis=uncertaintyAnalysis)
    # validate_argon(uncertaintyAnalysis=uncertaintyAnalysis)
    # validate_helium(uncertaintyAnalysis=uncertaintyAnalysis)

    # JT_helium_nitrogen(uncertaintyAnalysis=uncertaintyAnalysis)
    JT_helium_neon(uncertaintyAnalysis=uncertaintyAnalysis)

    plt.show()
