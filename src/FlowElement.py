from src.helper import *
import math

class FlowElement():
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



    def __repr__(self):
        string = "Flow Element at location x={:.2f} \n".format(self.x)
        for attr in ['M', 'p', 'rho', 'T']:
            val = self.__getattribute__(attr)
            string += '{} : {:.2f} \n'.format(attr, val)
        return string


    
def isentropic_propagate(fe : FlowElement, M1 = None, p1 = None, T1 = None, A1=None) -> FlowElement:
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
    M1 = fe.M
    p1 = fe.p
    T1 = fe.T
    gamma = fe.gamma

    # normal shock relations
    M2 = math.sqrt(((gamma - 1) * M1**2 + 2) / (2 * gamma * M1**2 - (gamma - 1))) if M1 > 1 else 1-1e-6
    p2 = p1 * (1 + (2 * gamma / (gamma + 1)) * (M1**2 - 1))
    T2 = T1 * ((p2 / p1) * ((gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2))

    post_shock_element = FlowElement(M2, p2, T2, gamma, fe.R, fe.A, fe.x)

    return post_shock_element
    





