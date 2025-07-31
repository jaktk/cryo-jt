import os
import sys
import time
import itertools
import csv
import numpy as np
import pandas as pd
import math
from scipy import optimize
from sklearn.metrics import r2_score
from numpy.polynomial.chebyshev import Chebyshev
import matplotlib as mpl
import matplotlib.pyplot as plt
import CoolProp
import CoolProp.CoolProp as CP
from FluidProps import FluidProps
from get_git_root import get_git_root

root = "~/codes/REFPROP"
if os.path.exists(root):
    CP.set_config_string(CoolProp.ALTERNATIVE_REFPROP_PATH, root)
else:
    print(f"{root} does not exist. Set path to REFPROP.")
    raise


class CernoxCal(object):
    serNum = ''

    def __init__(self, serNum):
        self.serNum = serNum
        self.coefficients = self.set_coefficients()
        self.tbl = self.set_tbl()
        self.dat = self.set_dat()

    def set_coefficients(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.cof')
        c_dict = {}
        with open(fpath, 'r') as fh:
            coefs = []
            limits = ''
            while True:
                line = fh.readline()
                if not line:
                    c_dict[limits] = {'Zlower': Zlower, 'Zupper': Zupper, 'coefs': coefs}
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
                    c_dict[limits] = {'Zlower': Zlower, 'Zupper': Zupper, 'coefs': coefs}
                    coefs = []
        return c_dict

    def set_tbl(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.tbl')
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def set_dat(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.dat')
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def get_T_from_coefs(self, R):
        T = 0
        for key, value in self.coefficients.items():
            L, U = key.split(',')
            L, U = float(L), float(U)
            if R >= L and R < U:
                k = ((np.log10(R) - value['Zlower']) - (value['Zupper'] - np.log10(R))) / (value['Zupper'] - value['Zlower'])
                for i, c in enumerate(value['coefs']):
                    T += c * np.cos(i * np.arccos(k))
        return T

    def get_T_from_tbl(self, R):
        f = interpolate.interp1d(self.tbl[('Resistance', '(Ohms)')], self.tbl[('Temp.', '(Kelvin)')])
        return f(R)

    def get_interp_abs_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Ti = self.get_T_from_tbl(R)
        return Tc - Ti

    def get_interp_rel_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Ti = self.get_T_from_tbl(R)
        return (Tc - Ti) / Ti

    def get_deltaT_dat_cof(self):
        T, DT = [], []
        for row in self.dat.iterrows():
            Tdat = row[1][('Temperature','(Kelvin)')]
            T.append(Tdat)
            DT.append(Tdat - self.get_T_from_coefs(row[1][('Resistance','(Ohms)')]))
        return T, DT

    def get_polynomial_fit_uncert(self, T):
        """ uncertainty in mK from fitting Chebyshev polynomial """
        d = {
        'X93303': [{'Tmin': 1.4 , 'Tmax': 14.1 , 'N': 31, 'n': 9, 'DTrms': 0.93},
                  {'Tmin': 14.1, 'Tmax': 80.0 , 'N': 31, 'n': 7, 'DTrms': 1.60},
                  {'Tmin': 80.0, 'Tmax': 325.0, 'N': 32, 'n': 8, 'DTrms': 5.84}],
        'X115143': [{'Tmin': 20.0, 'Tmax': 95.0, 'N': 28, 'n': 7, 'DTrms': 0.91},
                   {'Tmin': 95.0, 'Tmax': 325.0, 'N': 29, 'n': 9, 'DTrms': 3.45}],
        'X115888': [{'Tmin': 20.0, 'Tmax': 95.3, 'N': 28, 'n': 7, 'DTrms': 1.17},
                   {'Tmin': 95.3, 'Tmax': 325.0, 'N': 29, 'n': 8, 'DTrms': 4.06}]
        }
        for _ in d[self.serNum]:
            if T >= _['Tmin'] and T < _['Tmax']:
                sigm2 = _['N'] / (_['N'] - _['n']) * _['DTrms']**2
                break
            else:
                sigm2 = 0
        return 2 * sigm2**0.5


class TempUncertainty:
    sensor = ''
    def __init__(self, sensor):
        assert sensor in ['X93303','X115143','X115888']
        self.sensor = sensor
        self.cernox = CernoxCal(sensor)

    def get_R_from_T(self, T):
        dfx0 = self.cernox.dat.iloc[(self.cernox.dat[('Temperature','(Kelvin)')]-T).abs().argsort()[:1]]
        res = optimize.minimize(fun = lambda x: np.abs(self.cernox.get_T_from_coefs(x) - T),
                                x0 = dfx0[('Resistance','(Ohms)')])
        return res.x[0]

    def get_cernox_uncertainty(self, T):
        """ uncertainty in mK for temperature from 20 K to 300 K """
        return -5e-05 * T**2 + 0.1281 * T + 6.4603

    def get_cabtr_uncertainty(self, T):
        """ uncertainty in mK for temperature from 20 K to 300 K """
        return 2e-05 * T**2 + 0.0695 * T - 0.1001

    def get_polynomial_fit_uncert(self, T):
        """ uncertainty in mK from fitting Chebyshev polynomial """
        d = {
        'X93303': [{'Tmin': 1.4 , 'Tmax': 14.1 , 'N': 31, 'n': 9, 'DTrms': 0.93},
                  {'Tmin': 14.1, 'Tmax': 80.0 , 'N': 31, 'n': 7, 'DTrms': 1.60},
                  {'Tmin': 80.0, 'Tmax': 325.0, 'N': 32, 'n': 8, 'DTrms': 5.84}],
        'X115143': [{'Tmin': 20.0, 'Tmax': 95.0, 'N': 28, 'n': 7, 'DTrms': 0.91},
                   {'Tmin': 95.0, 'Tmax': 325.0, 'N': 29, 'n': 9, 'DTrms': 3.45}],
        'X115888': [{'Tmin': 20.0, 'Tmax': 95.3, 'N': 28, 'n': 7, 'DTrms': 1.17},
                   {'Tmin': 95.3, 'Tmax': 325.0, 'N': 29, 'n': 8, 'DTrms': 4.06}]
        }

        for i in d[self.sensor]:
            if T >= i['Tmin'] and T < i['Tmax']:
                sigm2 = i['N'] / (i['N'] - i['n']) * i['DTrms']**2
                break
            else:
                sigm2 = 0
        return 2 * sigm2**0.5

    def __call__(self, T):
        cabtr = self.get_cabtr_uncertainty(T)
        cernox = self.get_cernox_uncertainty(T)
        poly = self.get_polynomial_fit_uncert(T)
        return cabtr + cernox + poly


class Isenthalpic(object):
    def __init__(self, fluid):
        self.fluid = fluid
        self.AS = CP.AbstractState("HEOS", fluid)

    def get_hmass(self, p, T):
        self.AS.update(CP.PT_INPUTS, p*1e6, T)
        return self.AS.hmass()

    def get_T(self, p, h):
        self.AS.update(CP.HmassP_INPUTS, h, p*1e6)
        return self.AS.T()

    def get_JT_coefficient(self, p, T):
        """
        CoolProp does not provide a routine for calculating Joule-Thomson coefficient
        calculate it with derivatives from reduced Helmholtz energy
            return Joule-Thomson coefficient in K/MPa
        """
        self.AS.update(CP.PT_INPUTS, p*1e6, T)
        tau = self.AS.tau()
        delta = self.AS.delta()
        Tau2_d2alpha0_dTau2 = tau**2 * self.AS.d2alpha0_dTau2()
        Delta_dalphar_dDelta = delta * self.AS.dalphar_dDelta()
        Delta2_d2alphar_dDelta2 = delta**2 * self.AS.d2alphar_dDelta2()
        Tau2_d2alphar_dTau2 = tau**2 * self.AS.d2alphar_dTau2()
        Delta_Tau_d2alphar_dDelta_dTau = delta * tau * self.AS.d2alphar_dDelta_dTau()
        
        uJT_rhoR = - (Delta_dalphar_dDelta + Delta2_d2alphar_dDelta2 + Delta_Tau_d2alphar_dDelta_dTau) / (
            (1 + Delta_dalphar_dDelta - Delta_Tau_d2alphar_dDelta_dTau)**2  - (Tau2_d2alpha0_dTau2 +
                Tau2_d2alphar_dTau2) * (1 + 2 * Delta_dalphar_dDelta + Delta2_d2alphar_dDelta2))
        return uJT_rhoR / self.AS.rhomolar() / self.AS.gas_constant() * 1e6


def monte_carlo_pures(fluid, Tvec, pvec, fitterClass=Chebyshev, random_size=1000):
    FP = FluidProps(fluid)

    pvecmax, Tvecmax = max(pvec), Tvec[pvec.index(max(pvec))]
    p = np.arange(min(pvec), pvecmax+0.05, 0.05)
    hmass = FP.get_hmass(p=pvecmax, T=Tvecmax)
    T = np.array([FP.get_T(p=pi, h=hmass) for pi in p])
    Tvec = [FP.get_T(p=pi, h=hmass) for pi in pvec] # ideal values

    X93303 = TempUncertainty('X93303')
    dT = np.array([X93303(Ti) for Ti in Tvec]) * 1e-3 # K

    # find polynomial degree of the isenthalpic line
    for deg_h in range(0, 50):
        fit = fitterClass(fitterClass.fit(x=p, y=T, deg=deg_h).convert().coef)
        if r2_score(T, fit(p)) > 0.999:
            break

    # add random error (uniform distribution) to temperature and pressure measurements
    p_d, T_d = {}, {}
    for i in range(len(pvec)):
        p_d[i] = pvec[i] + np.random.uniform(low=-1.37e-3, high=1.37e-3, size=(random_size,))
        T_d[i] = Tvec[i] + np.random.uniform(low=-dT[i], high=+dT[i], size=(random_size,))

    # get pressure, temperature for measurements with randomized errors within sensors limit
    errSigm = {}
    for i in range(len(pvec)):
        errSigm[i] = []
    for j in range(random_size):
        prnd, Trnd = np.array([]), np.array([])
        for i in range(len(pvec)):
            prnd = np.append(prnd, p_d[i][j])
            Trnd = np.append(Trnd, T_d[i][j])

        # calculate JT coefficients and errors
        deg = len(pvec) - 2 if len(pvec) > 3 else len(pvec) - 1
        newSeries = fitterClass.fit(x=prnd, y=Trnd, deg=deg)
        fit = fitterClass(newSeries.convert().coef)
        newSeriesDeriv = fit.deriv(m=1)
        derivFit = fitterClass(newSeriesDeriv.convert().coef)

        uJT_meas = derivFit(prnd)
        uJT_calc = np.array([FP.get_JT_coefficient(p=_p, T=_T) for _p, _T in zip(prnd, Trnd)])
        _err = (uJT_meas - uJT_calc) / uJT_calc
        for i in range(len(pvec)):
            errSigm[i].append(_err[i])

    # calculate standard deviation for 95.45% confidence interval
    errStd = [2 * np.std(errSigm[i]) for i in range(len(pvec))]

    DTmax_Dpmax = (max(Tvec) - min(Tvec)) / (max(pvec) - min(pvec))
    return fluid, pvecmax, Tvecmax, DTmax_Dpmax, deg_h, len(pvec), errStd


def monte_carlo_mixtures(fluid, Tvec, pvec, x1=0.5, fitterClass=Chebyshev, random_size=1000, add_composition_uncert=True):
    FP = FluidProps(fluid)
    FP.set_composition_from_1st_fraction(x1)

    pvecmax, Tvecmax = max(pvec), Tvec[pvec.index(max(pvec))]
    p = np.arange(min(pvec), pvecmax+0.05, 0.05)
    hmass = FP.get_hmass(p=pvecmax, T=Tvecmax)
    T = np.array([FP.get_T(p=_, h=hmass) for _ in p])
    Tvec = [FP.get_T(p=_, h=hmass) for _ in pvec] # idealizing the measurements
    xvec = np.array([x1] * len(pvec))

    X93303 = TempUncertainty('X93303')
    dT = np.array([X93303(_) for _ in Tvec]) * 1e-3 # K

    # find polynomial degree of the h=const. line
    for deg_h in range(0, 50):
        fit = fitterClass(fitterClass.fit(x=p, y=T, deg=deg_h).convert().coef)
        if r2_score(T, fit(p)) > 0.999:
            break

    # add randomly distributed error to temperature, pressure, and composition
    p_d, T_d, x_d = {}, {}, {}
    for i in range(len(pvec)):
        p_d[i] = pvec[i] + np.random.uniform(low=-1.37e-3, high=1.37e-3, size=(random_size,))
        T_d[i] = Tvec[i] + np.random.uniform(low=-dT[i], high=+dT[i], size=(random_size,))
        if add_composition_uncert:
            if x1 <= 1e-3:
                limLow, limHigh = -x1, 1e-3
            elif x1 >= 0.999:
                limLow, limHigh = -1e-3, (1-x1)
            else:
                limLow, limHigh = -1e-3, 1e-3
            x_d[i] = xvec[i] + np.random.uniform(low=limLow, high=limHigh, size=(random_size,))
        else:
            x_d[i] = [xvec[i]] * random_size

    # get pressure, temperature, composition with randomized error
    errSigm = {}
    for i in range(len(pvec)):
        errSigm[i] = []
    for j in range(random_size):
        prnd, Trnd, xrnd = np.array([]), np.array([]), np.array([])
        for i in range(len(pvec)):
            prnd = np.append(prnd, p_d[i][j])
            Trnd = np.append(Trnd, T_d[i][j])
            xrnd = np.append(xrnd, x_d[i][j])

        # calculate JT coefficients and errors
        deg = len(pvec) - 2 if len(pvec) > 3 else len(pvec) - 1
        newSeries = fitterClass.fit(x=prnd, y=Trnd, deg=deg)
        fit = fitterClass(newSeries.convert().coef)
        newSeriesDeriv = fit.deriv(m=1)
        derivFit = fitterClass(newSeriesDeriv.convert().coef)

        uJT_meas = derivFit(prnd)
        uJT_calc = np.array([])
        for _p, _T, _x in zip(prnd, Trnd, xrnd):
            FP.set_composition_from_1st_fraction(_x)
            uJT_calc = np.append(uJT_calc, FP.get_JT_coeffcient(p=_p, T=_T))
        _err = (uJT_meas - uJT_calc) / uJT_calc

        if (-1 in np.around(_err, 3)) or (1 in np.around(_err, 3)):
            _err = np.array([0] * len(_err))

        for i in range(len(pvec)):
            errSigm[i].append(_err[i])

    # calculate standard deviation for 95.45% confidence interval
    errStd = [2 * np.std(errSigm[i]) for i in range(len(pvec))]

    DTmax_Dpmax = (max(Tvec) - min(Tvec)) / (max(pvec) - min(pvec))
    return fluid, pvecmax, Tvecmax, x1, DTmax_Dpmax, deg_h, len(pvec), errStd


def monte_carlo_pures_general(fluid, T1, p1, npts, p2=0.1, fitterClass=Chebyshev, random_size=1000):
    """
    routine for calculating the uncertainty for pure fluids
    from p1, T1, p2 (floats) and npts (int)
    """
    dp = p1 - p2
    pvec = list(np.linspace(p2 + 0.1*dp, p2 + 0.9*dp, num=npts))

    FP = FluidProps(fluid)
    hmass = FP.get_hmass(p=p1, T=T1)
    Tvec = np.array([FP.get_T(p=_, h=hmass) for _ in pvec])

    res = monte_carlo_pures(fluid=fluid, Tvec=Tvec, pvec=pvec, fitterClass=fitterClass, random_size=random_size)

    fluid, p2max, T2max, DTmax_DPmax, deg_h, npts, errStd = res
    return fluid, p1, T1, DTmax_DPmax, deg_h, npts, errStd


res = monte_carlo_pures_general('Nitrogen', 160, 8, 5)
print(np.array(res[-1]), np.mean(np.array(res[-1])[1:-1]))

res = monte_carlo_pures_general('Nitrogen', 160, 8, 6)
print(np.array(res[-1]), np.mean(np.array(res[-1])[1:-1]))


def monte_carlo_mixtures_general(fluid, T1, p1, npts, p2=0.1, x1=0.5, fitterClass=Chebyshev, random_size=1000, add_composition_uncert=True):
    """
    routine for calculating the uncertainty for mixtures
    from p1, T1, p2 (floats) and npts (int)
    """
    dp = p1 - p2
    pvec = np.linspace(p2 + 0.1*dp, p2 + 0.9*dp, num=npts)

    FP = FluidProps(fluid)
    hmass = FP.get_hmass(p=p1, T=T1)
    Tvec = np.array([FP.get_T(p=_, h=hmass) for _ in pvec])

    res = monte_carlo_mixtures(fluid=fluid, Tvec=Tvec, pvec=pvec, x1=x1, fitterClass=Chebyshev,
        random_size=random_size, add_composition_uncert=add_composition_uncert)

    fluid, pvecmax, Tvecmax, x1, DTmax_Dpmax, deg_h, npts, errStd = res
    return fluid, p1, T1, x1, DTmax_Dpmax, deg_h, npts, errStd


def time_string(t):
    if t < 60:
        return t, 's'
    elif t >= 60 and t < 3600:
        return t/60, 'min'
    else:
        return t/3600, 'h'


def dump_monte_carlo_pures_overview(fileNum, fluid, Tarray, parray, narray, fitterClass):
    print('Total number of iterations: {}'.format(len(Tarray) * len(parray) * len(narray)))
    firstTic = time.time()
    res = {}
    i = 0
    for T in Tarray:
        for p in parray:
            for npts in narray:
                try:
                    tic = time.time()
                    res[i] = monte_carlo_pures_general(fluid = fluid,
                                                       p1 = p,
                                                       T1 = T,
                                                       npts = npts,
                                                       random_size = 200 * npts,
                                                       fitterClass = fitterClass)
                    toc = time.time()
                    t, u = time_string(toc-tic)
                    ttot, utot = time_string(toc-firstTic)
                    print('iter {}: {}, finished in {:.2f} {}. Total time elapsed: {:.2f} {}'.format(i, fluid, t, u, ttot, utot))
                except:
                    print('iter {}: {} FAILED. T = {} K, p = {} MPa, n = {}'.format(i, fluid, T, p, npts))
                i += 1
    df = pd.DataFrame.from_dict(res,
        orient='index',
        columns=['fluid','p1/MPa','T1/K','DT/Dp/(K/MPa)','deg_h=const','n_pts','ErrArr'])
    _ = str(fitterClass).rstrip('\'>').split('.')[-1]
    df.to_csv('slope_error_summary_monte_carlo_{}_{}_{}.csv'.format(fluid, _, fileNum), sep=',', index=False)


def dump_monte_carlo_mixtures_overview(fileNum, fluid, xarray, Tarray, parray, narray):
    print('Total number of iterations: {}'.format(len(xarray) * len(Tarray) * len(parray) * len(narray)))
    firstTic = time.time()
    res = {}
    i = 0
    for x in xarray:
        for T in Tarray:
            for p in parray:
                for npts in narray:
                    try:
                        tic = time.time()
                        resTrue = monte_carlo_mixtures_general(fluid = fluid,
                                                               x1 = x,
                                                               p1 = p,
                                                               T1 = T,
                                                               npts = npts,
                                                               random_size = 200 * npts,
                                                               add_composition_uncert = True)
                        resFalse = monte_carlo_mixtures_general(fluid = fluid,
                                                                x1 = x,
                                                                p1 = p,
                                                                T1 = T,
                                                                npts = npts,
                                                                random_size = 200 * npts,
                                                                add_composition_uncert = False)
                        res[i] = {'fluid': resTrue[0],
                                  'p1/MPa': resTrue[1],
                                  'T1/K': resTrue[2],
                                  'x1': resTrue[3],
                                  'DT/Dp/(K/MPa)': resTrue[4],
                                  'deg_h=const': resTrue[5],
                                  'n_pts': resTrue[6],
                                  'ErrArrWithCmpUncert': resTrue[7],
                                  'ErrArrNoCmpUncert': resFalse[7]}
                        toc = time.time()
                        t, u = time_string(toc-tic)
                        ttot, utot = time_string(toc-firstTic)
                        print('iter {}: {}, finished in {:.2f} {}. Total time elapsed: {:.2f} {}'.format(i, fluid, t, u, ttot, utot))
                    except ValueError:
                        print('iter {}: {} FAILED. T = {} K, p = {} MPa, n = {}'.format(i, fluid, T, p, npts))
                    i += 1
    df = pd.DataFrame.from_dict(res,
        orient='index',
        columns=['fluid','p1/MPa','T1/K','x1','DT/Dp/(K/MPa)','deg_h=const','n_pts','ErrArrWithCmpUncert','ErrArrNoCmpUncert'])
    df.to_csv('slope_error_summary_monte_carlo_{}_{}.csv'.format(fluid[0]+'-'+fluid[1], fileNum), sep=',', index=False)


def file_post_process(filename):
    with open(filename, mode='r') as f:
        lines = [line.rstrip().replace('"','') for line in f]
    mem = ''
    _filename = filename.replace('.csv','')+'_POST.csv'
    with open(_filename, mode='w+') as f:
        f.write(lines[0]+'\n')
        for line in lines[1:]:
            if line[-1] == ']':
                if mem == '':
                    f.write(line+'\n')
                else:
                    f.write(mem+line+'\n')
                    mem = ''
            else:
                mem += line


def monte_carlo_N2():
    files = ["Nitrogen_150K_6MPa.csv", "Nitrogen_160K_12MPa.csv",
             "Nitrogen_160K_5MPa.csv", "Nitrogen_160K_6MPa.csv",
             "Nitrogen_160K_7MPa.csv", "Nitrogen_160K_9MPa.csv"]
    measurements = []
    for file in files:
        df = os.path.join(get_git_root(os.getcwd()), "data", "p_T_pairs", file)
        measurements.append({"pin": df["PT101/MPa"].mean(axis=0),
                             "Tin": df["TT101/K"].mean(axis=0),
                             "pout": list(df["PT102/MPa"])
                             "Tout": list(df["TT102/K"])})

    for m in measurements:
        while True:
            try:
                res = monte_carlo_pures(fluid = 'Nitrogen',
                                        Tvec = m['Tout'],
                                        pvec = m['pout'],
                                        fitterClass = Chebyshev,
                                        random_size = 10000)
                break
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                continue
        print('pin: {}, Tin: {}, deg_h: {}, err: {}'.format(m['pin'], m['Tin'], res[0], res[1][1:-1]))


def monte_carlo_Ar():
    measurements = [
        {
            'pin': 12.0,
            'Tin': 180,
            'pout': [11.01050034, 9.021099854, 8.045300293, 7.020600128, 5.979700089, 5.014699936],
            'Tout': [177.8821259, 171.7108459, 168.0425415, 163.5248413, 158.1468811, 152.1369629]},
        {
            'pin': 5.0,
            'Tin': 180.0,
            'pout': [3.990000153, 2.998999977, 1.973800087, 1.007999992],
            'Tout': [174.2562714, 165.6718597, 155.6197662, 144.6188965]}
    ]

    for m in measurements:
        while True:
            try:
                res = monte_carlo_pures(fluid = 'Argon',
                                        Tvec = m['Tout'],
                                        pvec = m['pout'],
                                        fitterClass = Chebyshev,
                                        random_size = 10000)
                break
            except:
                print('failed')
                continue
        print('pin: {}, Tin: {}, deg_h: {}, err: {}'.format(m['pin'], m['Tin'], res[0], res[1][1:-1]))


def monte_carlo_He():
    measurements = [
        {
            'pin': 7.0,
            'Tin': 65.0,
            'pout': [5.004700089, 3.918799973, 3.083499908, 2.018099976],
            'Tout': [67.40324097, 67.63872986, 68.28799744, 68.42211517]
        },
        {
            'pin': 5.0,
            'Tin': 64.99,
            'pout': [4.042599869, 2.99829998,  2.001300049, 1.004399967],
            'Tout': [66.38734436, 66.78726654, 67.05940704, 67.38295135]
        },
        {
            'pin': 10.0,
            'Tin': 140.6,
            'pout': [9.022100067, 8.023799896, 5.975899887, 3.179000092],
            'Tout': [141.8423767, 142.219635,  143.2433319, 144.7186127]
        }
    ]

    for m in measurements:
        while True:
            try:
                res = monte_carlo_pures(fluid = 'Helium',
                                        Tvec = m['Tout'],
                                        pvec = m['pout'],
                                        fitterClass = Chebyshev,
                                        random_size = 10000)
                break
            except:
                print('failed')
                continue
        print('pin: {}, Tin: {}, deg_h: {}, err: {}'.format(m['pin'], m['Tin'], res[0], res[1][1:-1]))


def monte_carlo_HeNe():
    measurements = [
        {
            'pin': 10.0,
            'Tin': 65.70153,
            'x1': 0.4623,
            'pout': [7.998799896, 6.073199844, 4.154499817, 2.204599953],
            'Tout': [65.48348999, 64.09143829, 62.50511551, 60.45336533]
        },
        {
            'pin': 10.0002,
            'Tin': 79.899048,
            'x1': 0.4666,
            'pout': [6.526999664, 5.99090004,  5.005799866, 4.005099869],
            'Tout': [79.49318695, 79.15239655, 78.55823517, 78.0082016]
        },
        {
            'pin': 7.9997,
            'Tin': 64.99115,
            'x1': 0.4673,
            'pout': [6.065599823, 5.009549968, 3.87179985,  2.977099991],
            'Tout': [64.5326004,  63.6618309,  62.62127304, 61.67658615]
        },
        {
            'pin': 5.0,
            'Tin': 80.00457,
            'x1': 0.467, 
            'pout': [3.909799957, 3.024099922, 1.943799973, 0.994699955],
            'Tout': [80.71064758, 79.9641037,  79.1587677,  78.42686554]
        },
        {
            'pin': 7.0003,
            'Tin': 79.994446,
            'x1': 0.4667,
            'pout': [5.018399811, 4.001399994, 3.004999924],
            'Tout': [79.87950134, 79.25732422, 78.63114166]
        },
        {
            'pin': 10.0001,
            'Tin': 65.142357,
            'x1': 0.3138,
            'pout': [9.07460022, 8.152199936, 5.967166742, 3.975699997, 2.037100029],
            'Tout': [65.41204071, 64.4666748, 62.02532959, 59.168396,   55.79568863]
        },
        {
            'pin': 7.9999,
            'Tin': 80.435432,
            'x1': 0.3264,
            'pout': [6.974299622, 6.036999893, 3.975099945, 1.992499924],
            'Tout': [80.75549316, 79.88356781, 78.01207733, 75.90357208]
        },
        {
            'pin': 5.0006,
            'Tin': 65.161682,
            'x1': 0.3279,
            'pout': [3.134933344, 2.091350079, 1.047200012],
            'Tout': [63.78694153, 61.97766876, 59.99108887]
        },
        {
            'pin': 5.0,
            'Tin': 65.013657,
            'x1': 0.3967,
            'pout': [3.843500137, 3.135400009, 2.072699928, 0.980300045],
            'Tout': [65.14979553, 64.08061981, 62.58144379, 60.87456131]
        },
        {
            'pin': 7.0085,
            'Tin': 64.994644,
            'x1': 0.3942,
            'pout': [5.993350029, 4.964900017, 4.055400085, 3.070700073, 1.999699974],
            'Tout': [65.26316071, 64.01760101, 62.85807037, 61.49729919, 60.02370453]
        },
        {
            'pin': 5.0,
            'Tin': 65.052849,
            'x1': 0.2154,
            'pout': [3.954000092, 3.109225035, 2.044149971, 1.033400011, 0.146700001],
            'Tout': [65.03496552, 63.12421417, 60.68468857, 58.11508179, 55.57246399]
        },
        {
            'pin': 7.0,
            'Tin': 65.047478,
            'x1': 0.2005,
            'pout': [6.011100006, 5.386174965, 4.028900146, 3.272700119, 1.984300041, 1.003100014, 0.175800002],
            'Tout': [65.1322403,  63.9445076,  61.42861938, 59.8692894,  56.83626175, 54.13127899, 51.43774033]
        }
    ]

    for m in measurements:
        while True:
            try:
                res = monte_carlo_mixtures(fluid = ('Helium', 'Neon'),
                                           Tvec = m['Tout'],
                                           pvec = m['pout'],
                                           x1 = m['x1'],
                                           random_size = 10000,
                                           add_composition_uncert = True)
                break
            except:
                print('failed')
                continue
        print('pin: {}, Tin: {}, x1: {}, deg_h: {}, err: {}'.format(m['pin'], m['Tin'], m['x1'], res[0], res[1][1:-1]))


def monte_carlo_HeN2():
    measurements = [
        {
            'pin': 7.9999,
            'Tin': 160.008926,
            'x1': 0.1451,
            'pout': [5.99640007,  5.05929985,  2.991799927, 1.943700027, 0.993999958],
            'Tout': [153.3766632, 149.3185425, 138.3912659, 131.4246979, 124.1477966]
        },
        {
            'pin': 5.0,
            'Tin': 160.000748,
            'x1': 0.1449,
            'pout': [3.980599976, 2.013899994, 1.000399971],
            'Tout': [156.9343109, 146.4505005, 140.2199402]
        },
        {
            'pin': 5.0,
            'Tin': 139.999939,
            'x1': 0.145,
            'pout': [3.997700119, 3.006200027, 2.005900002, 1.039000034],
            'Tout': [135.7933197, 129.3087616, 121.5911865, 112.7962875]
        },
        {
            'pin': 5.0,
            'Tin': 140.005066,
            'x1': 0.5025,
            'pout': [3.998400116, 3.276900101, 2.009799957, 0.990100002],
            'Tout': [139.6854858, 137.325943,  134.9553375, 132.3852539]
        },
        {
            'pin': 8.0,
            'Tin': 140.00708,
            'x1': 0.5026,
            'pout': [6.950800323, 5.979199982, 5.00890007,  4.007300186, 2.018499947, 1.00340004],
            'Tout': [139.530899,  137.6961823, 135.8578796, 133.8177643, 129.152832,  126.3856506]
        }
    ]

    for m in measurements:
        while True:
            try:
                res = monte_carlo_mixtures(fluid = ('Helium', 'Nitrogen'),
                                           Tvec = m['Tout'],
                                           pvec = m['pout'],
                                           x1 = m['x1'],
                                           random_size = 10000,
                                           add_composition_uncert = True)
                break
            except:
                print('failed')
                continue
        print('pin: {}, Tin: {}, x1: {}, deg_h: {}, err: {}'.format(m['pin'], m['Tin'], m['x1'], res[0], res[1][1:-1]))


if __name__ == '__main__':
    # Chebyshev, Polynomial, Legendre, Laguerre, Hermite
    fitterClass = Chebyshev
    fileNum = 0
    # dropExtremities = False

    # file_post_process('slope_error_summary_monte_carlo_Argon_Chebyshev_0.csv')
    # file_post_process('slope_error_summary_monte_carlo_Helium_Chebyshev_0.csv')
    # file_post_process('slope_error_summary_monte_carlo_Neon_Chebyshev_0.csv')
    # file_post_process('slope_error_summary_monte_carlo_Nitrogen_Chebyshev_0.csv')
    
    """ Argon: 576 iterations """
    # dump_error_summary(fileNum = fileNum,
    #                    fluid = 'Argon',
    #                    Tarray = np.arange(195, 255, 5),
    #                    parray = np.arange(5, 13, 1),
    #                    narray = np.arange(3, 9, 1),
    #                    fitterClass = fitterClass,
    #                    dropExtremities = dropExtremities)
    # dump_monte_carlo_pures_overview(fileNum = fileNum,
    #                                 fluid = 'Argon',
    #                                 Tarray = np.arange(195, 255, 5),
    #                                 parray = np.arange(5, 13, 1),
    #                                 narray = np.arange(3, 9, 1),
    #                                 fitterClass = fitterClass)

    """ Nitrogen: 480 iterations """
    # dump_error_summary(fileNum = fileNum,
    #                    fluid = 'Nitrogen',
    #                    Tarray = np.arange(160, 260, 10),
    #                    parray = np.arange(5, 13, 1),
    #                    narray = np.arange(3, 9, 1),
    #                    fitterClass = fitterClass,
    #                    dropExtremities = dropExtremities)
    # dump_monte_carlo_pures_overview(fileNum = fileNum,
    #                                 fluid = 'Nitrogen',
    #                                 Tarray = np.arange(160, 260, 10),
    #                                 parray = np.arange(5, 13, 1),
    #                                 narray = np.arange(3, 9, 1),
    #                                 fitterClass = fitterClass)

    """ Helium: 180 iterations """
    # dump_error_summary(fileNum = fileNum,
    #                    fluid = 'Helium',
    #                    Tarray = np.arange(60, 160, 10),
    #                    parray = np.array([5, 8, 12]),
    #                    narray = np.arange(3, 9, 1),
    #                    fitterClass = fitterClass,
    #                    dropExtremities = dropExtremities)
    # dump_monte_carlo_pures_overview(fileNum = fileNum,
    #                                 fluid = 'Helium',
    #                                 Tarray = np.arange(60, 160, 10),
    #                                 parray = np.array([5, 8, 12]),
    #                                 narray = np.arange(3, 9, 1),
    #                                 fitterClass = fitterClass)

    """ Neon: 576 iterations """
    # dump_error_summary(fileNum = fileNum,
    #                    fluid = 'Neon',
    #                    Tarray = np.arange(65, 125, 5),
    #                    parray = np.arange(5, 13, 1),
    #                    narray = np.arange(3, 9, 1),
    #                    fitterClass = fitterClass,
    #                    dropExtremities = dropExtremities)
    # dump_monte_carlo_pures_overview(fileNum = fileNum,
    #                                 fluid = 'Neon',
    #                                 Tarray = np.arange(65, 125, 5),
    #                                 parray = np.arange(5, 13, 1),
    #                                 narray = np.arange(3, 9, 1),
    #                                 fitterClass = fitterClass)


    """ Helium-Neon: 5760 iterations """
    # dump_error_summary_mixture(fileNum = fileNum,
    #                            fluid = ('Helium','Neon'),
    #                            xarray = np.arange(0, 1.1, 0.1),
    #                            Tarray = np.arange(65, 125, 5),
    #                            parray = np.arange(5, 13, 1),
    #                            narray = np.arange(3, 9, 1))
    # dump_monte_carlo_mixtures(fileNum = fileNum,
    #                           fluid = ('Helium','Neon'),
    #                           xarray = np.arange(0, 1.1, 0.1),
    #                           Tarray = np.arange(65, 125, 5),
    #                           parray = np.arange(5, 13, 1),
    #                           narray = np.arange(3, 9, 1))

    # dump_monte_carlo_mixtures(fileNum = 0,
    #                           fluid = ('Helium','Neon'),
    #                           xarray = np.array([0.0, 0.1, 0.2]),
    #                           Tarray = np.arange(65, 125, 5),
    #                           parray = np.arange(5, 13, 1),
    #                           narray = np.arange(3, 9, 1))
    # dump_monte_carlo_mixtures(fileNum = 1,
    #                           fluid = ('Helium','Neon'),
    #                           xarray = np.array([0.3, 0.4, 0.5]),
    #                           Tarray = np.arange(65, 125, 5),
    #                           parray = np.arange(5, 13, 1),
    #                           narray = np.arange(3, 9, 1))
    # dump_monte_carlo_mixtures(fileNum = 2,
    #                           fluid = ('Helium','Neon'),
    #                           xarray = np.array([0.6, 0.7, 0.8]),
    #                           Tarray = np.arange(65, 125, 5),
    #                           parray = np.arange(5, 13, 1),
    #                           narray = np.arange(3, 9, 1))
    # dump_monte_carlo_mixtures(fileNum = 3,
    #                           fluid = ('Helium','Neon'),
    #                           xarray = np.array([0.9, 1.0]),
    #                           Tarray = np.arange(65, 125, 5),
    #                           parray = np.arange(5, 13, 1),
    #                           narray = np.arange(3, 9, 1))

    # uJT_error_mixture_simple(fluid=('Helium','Neon'), T1=80, p1=10, x1=0.8)
    # uJT_error_mixture_simple(fluid='Nitrogen', T1=160, p1=10)

    # monte_carlo_N2()
    # monte_carlo_Ar()
    # monte_carlo_He()
    # monte_carlo_HeNe()
    # monte_carlo_HeN2()

    plt.show()