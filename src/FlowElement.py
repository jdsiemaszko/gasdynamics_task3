"""
Useful class and function definitions for Task 3
"""
from src.helper import *
import math

class FlowElement():
    """
    Main computational element - a point in the flow.
    Can be propagated along isentropes, Fanno lines, and Rayleigh lines
    state variables: M, p, T
    other variables not stored, but can be computed as properties
    """
    def __init__(self, M, p, T, gamma=1.4, R=287.05, A = None, x = 0.0):
        
        self.__M = M
        self.__p = p
        self.__T = T
        self.__gamma = gamma
        self.__R = R
        self.__A = A
        self.__x = x

    @property
    def x(self):
        return self.__x

    @property
    def A(self):
        return self.__A

    @property
    def Dh(self):
        if self.A is None:
            return None
        return math.sqrt(self.A * 4 / math.pi)

    @property
    def gamma(self):
        return self.__gamma
    
    @property
    def R(self):
        return self.__R

    @property
    def p(self):
        return self.__p

    @property
    def M(self):
        return self.__M

    @property
    def T(self):
        return self.__T

    @property
    def rho(self):
        return self.p / self.T / self.R
    
    @property
    def a(self):
        return math.sqrt(self.gamma * self.R * self.T)
    
    @property
    def V(self):
        return self.a * self.M

    @property
    def Ttot(self):
        return self.T * (1 + ((self.gamma - 1) / 2) * self.M**2)

    @property
    def ptot(self):
        return self.p * (1 + ((self.gamma - 1) / 2) * self.M**2) ** (self.gamma / (self.gamma - 1))

    @property
    def rhotot(self):
        return self.rho * (1 + ((self.gamma - 1) / 2) * self.M**2) ** (1 / (self.gamma - 1))

    @property
    def M_max_T(self):
        return 1 / math.sqrt(self.gamma)

    @property
    def Cp(self):
        return self.R * self.gamma / (self.gamma-1)

    def __repr__(self):
        # helper for printing results (can be adjusted to print total properties, etc.)

        string = "Flow Element at location x={:.3f} \n".format(self.x)
        for attr in ['M', 'p', 'rho', 'T']:
            val = self.__getattribute__(attr)
            string += '{} : {:.2f} \n'.format(attr, val)
        return string


    
def isentropic_propagate(fe : FlowElement, M1 = None, p1 = None, T1 = None, A1=None) -> FlowElement:
    """
    Propagate FlowElement fe along an isentrope
    :param fe: a FlowElement instance
    :param M1: resulting Mach number, used to compute properties of the new point
    :param p1: resulting pressure (not implemented)
    :param T1: resulting temperature (not implemented)
    :param A1: resulting duct area, used to compute properties of the new point
    :return: fe1 -  the resulting FlowElement instance
    """
    if M1 is not None:
        p1 = isentropic_pressure_ratio(M1, fe.M, fe.gamma) * fe.p
        T1 = isentropic_temp_ratio(M1, fe.M, fe.gamma) * fe.T

        if fe.M == 0:
            A1 = fe.A
        else:
            A1 = isentropic_area_ratio(M1, fe.M, fe.gamma) * fe.A

        fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, A1, fe.x)

        return fe1
    elif A1 is not None:
        M1 = isentropic_Mach_from_area_ratio(A1 / fe.A, fe.M, fe.gamma)
        p1 = isentropic_pressure_ratio(M1, fe.M, fe.gamma) * fe.p
        T1 = isentropic_temp_ratio(M1, fe.M, fe.gamma) * fe.T
        fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, A1, fe.x)

        return fe1
    else:
        raise NotImplementedError('bad input')


def Fanno_propagate(fe: FlowElement, friction, distance=None, M1=None) -> FlowElement:
    """
        Propagate FlowElement fe along a Fanno line
        :param fe: a FlowElement instance
        :param M1: resulting Mach number, used to compute properties of the new point
        :param distance: length of the duct
        :return: fe1 -  the resulting FlowElement instance
        """
    if distance is not None:

        xfinal = distance + fe.x

        M1, _, _ = Fanno_integrator(fe.M, fe.p, fe.T, xfinal, fe.Dh, fe.gamma, friction, x_initial=fe.x)
        # M1 = Fanno_Mach_integrator(fe.M, xfinal, fe.Dh, fe.gamma, friction, x_initial=fe.x)
        p1 = Fanno_pressure_ratio(M1, fe.M, fe.gamma) * fe.p
        T1 = Fanno_temperature_ratio(M1, fe.M, fe.gamma) * fe.T

        fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, fe.A, fe.x + distance)

        return fe1

    elif M1 is not None:
        dist = Fanno_length_from_Mach_difference(M0=fe.M, Mfinal=M1, gamma=fe.gamma, friction=friction, Dh=fe.Dh)
        p1 = Fanno_pressure_ratio(M1, fe.M, fe.gamma) * fe.p
        T1 = Fanno_temperature_ratio(M1, fe.M, fe.gamma) * fe.T

        fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, fe.A, fe.x + dist)

        return fe1
    else:
        raise NotImplementedError('bad input')

def normal_shock(fe: FlowElement) -> FlowElement:
    """
    Normal shock relation
    :param fe: upstream FlowElement (should be supersonic)
    :return: downstream FlowElement (always subsonic)
    """
    M1 = fe.M
    if M1 < 1.0:
        return fe # no shock!

    p1 = fe.p
    T1 = fe.T
    gamma = fe.gamma

    # normal shock relations
    M2 = math.sqrt(((gamma - 1) * M1**2 + 2) / (2 * gamma * M1**2 - (gamma - 1))) if M1 > 1 else 1-1e-6
    p2 = p1 * (1 + (2 * gamma / (gamma + 1)) * (M1**2 - 1))
    T2 = T1 * ((p2 / p1) * ((gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2))

    post_shock_element = FlowElement(M2, p2, T2, gamma, fe.R, fe.A, fe.x)

    return post_shock_element

def Rayleigh_heat_addition(fp:FlowElement, dq)->FlowElement:
    """
    propagate along a Rayleigh line due to heat addition dq
    :param fp: initial FlowElement instance
    :param dq: heat addition (positive-> heat added to the flow)
    :return: final FlowElement instance
    """

    delta_Ttot = dq / fp.Cp
    Ttot2 = fp.Ttot + delta_Ttot

    M2 = Rayleigh_Mach_from_Ttot_ratio(Ttot2, fp.Ttot, fp.M, fp.gamma)
    p2 = Rayleigh_p_ratio(M2, fp.M, fp.gamma) * fp.p
    T2 = Rayleigh_T_ratio(M2, fp.M, fp.gamma) * fp.T

    return FlowElement(M2, p2, T2, A = fp.A)

def Rayleigh_propagate(fp, M2):
    """
        propagate along a Rayleigh line up to Mach number M2
        :param fp: initial FlowElement instance
        :param M2: final Mach number
        :return: final FlowElement instance
    """
    p2 = Rayleigh_p_ratio(M2, fp.M, fp.gamma) * fp.p
    T2 = Rayleigh_T_ratio(M2, fp.M, fp.gamma) * fp.T

    return FlowElement(M2, p2, T2, A=fp.A)








