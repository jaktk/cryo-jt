import sys
import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary


class FluidProps(object):
    """
    small REFPROP wrapper to simplify syntax in the main code
    """
    def __init__(self, fluid):
        self.set_root()
        self.RP = REFPROPFunctionLibrary(self.root)
        self.RP.SETPATHdll(self.root)
        self.MASS_BASE_SI = self.RP.GETENUMdll(0, "MASS BASE SI").iEnum
        self.MOLAR_BASE_SI = self.RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
        self.fluid = self.set_fluid(fluid)
        self.z = {1.0}

    def set_root(self):
        # Prefer REFPROP's standard RPPREFIX environment variable; fall back to
        # the conventional per-platform install locations.
        env_root = os.environ.get('RPPREFIX')
        candidates = [env_root] if env_root else []
        home_dir = os.environ.get('HOME', '')
        platform_default = {
            'linux': '/etc/REFPROP',
            'win32': 'c:/Program Files (x86)/REFPROP',
            'darwin': os.path.join(home_dir, 'codes', 'REFPROP'),
        }.get(sys.platform)
        if platform_default:
            candidates.append(platform_default)

        for root in candidates:
            if root and os.path.exists(root):
                self.root = root
                os.environ['RPPREFIX'] = root
                return

        raise FileNotFoundError(
            "REFPROP installation not found. Set the RPPREFIX environment "
            "variable to your REFPROP directory (tried: "
            f"{', '.join(c for c in candidates if c) or 'no candidates'})."
        )

    def set_fluid(self, fluid):
        return fluid if type(fluid) == str else '{};{}'.format(*fluid)

    def set_composition(self, z):
        self.z = z

    def set_composition_from_1st_fraction(self, x):
        """
        set mixture composition from the first component's concentration
        """
        self.z = [round(x, 8), round(1-x, 8)]

    def get_compressibility_factor(self, p, T):
        """
        return compressibility factor (dimensionless) as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'Z',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_p(self, rho, T):
        """
        return pressure in MPa as a function of:
            `rho` density in kg/m^3
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'DT', 'P',
                                 self.MASS_BASE_SI, 0, 0,
                                 rho, T, self.z)
        return res.Output[0] / 1e6

    def get_rhomass(self, p, T):
        """
        return density in kg/m^3 as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'D',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_rhomolar(self, p, T):
        """
        return density in mol/dm^3 as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'D',
                                 self.MOLAR_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_cp(self, p, T):
        """
        return specific heat at constant pressure in J/kg/K as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'CP',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_specific_heat_ratio(self, p, T):
        """
        return the specific heat ratio as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'CP/CV',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_hmass(self, p, T):
        """
        return enthalpy in J/kg as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'H',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_T(self, p, h):
        """
        return temperature in K as a function of:
            `p` pressure in MPa
            `h` enthalpy in J/kg
        """
        res = self.RP.REFPROPdll(self.fluid, 'PH', 'T',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, h, self.z)
        return res.Output[0]

    def get_critical_T(self):
        """
        return critical temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PQ', 'TC',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0]

    def get_critical_p(self):
        """
        return critical pressure in MPa
        """
        res = self.RP.REFPROPdll(self.fluid, 'TQ', 'PC',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0] / 1e6

    def get_triple_T(self):
        """
        return triple point temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PQ', 'TTRP',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0]

    def get_triple_p(self):
        """
        return triple point pressure in MPa
        """
        res = self.RP.REFPROPdll(self.fluid, 'TQ', 'PTRP',
                                 self.MASS_BASE_SI, 0, 0,
                                 1, 1, self.z)
        return res.Output[0] / 1e6

    def get_saturation_T(self, p):
        """
        return saturation temperature in K as a function of:
            `p` pressure in MPa
        """
        res = self.RP.REFPROPdll(self.fluid, 'PQ', 'T',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, 1, self.z)
        return res.Output[0]

    def get_saturation_p(self, T):
        """
        return saturation pressure in MPa as a function of:
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'TQ', 'P',
                                 self.MASS_BASE_SI, 0, 0,
                                 T, 1, self.z)
        return res.Output[0] / 1e6

    def get_speed_of_sound(self, p, T):
        """
        return speed of sound in m/s as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'W',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0]

    def get_JT_coefficient(self, p, T):
        """
        return Joule-Thomson coefficient in K/MPa as a function of:
            `p` pressure in MPa
            `T` temperature in K
        """
        res = self.RP.REFPROPdll(self.fluid, 'PT', 'JT',
                                 self.MASS_BASE_SI, 0, 0,
                                 p*1e6, T, self.z)
        return res.Output[0] * 1e6


def test():
    """
    test the wrapper with helium-neon EOS
    """
    FP = FluidProps(('Helium', 'Neon'))
    FP.set_composition([0.5, 0.5])
    dh = FP.get_hmass(p=0.1, T=300) - FP.get_hmass(p=0.1, T=100)
    print('dh = {:.5f} kJ/kg'.format(dh/1e3))
    print('u_JT = {:.5f} K/MPa'.format(FP.get_JT_coefficient(p=1, T=60)))


if __name__ == '__main__':
    test()
