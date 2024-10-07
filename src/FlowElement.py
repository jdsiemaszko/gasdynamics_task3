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
    
def isentropic_propagate(fe : FlowElement, M1 = None, p1 = None, T1 = None) -> FlowElement:
    if M1 is not None:
        p1 = isentropic_p_over_p_tot(fe.M, fe.gamma)
        T1 = isentropic_T_over_T_tot(fe.M, fe.gamma)

        fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, fe.A, fe.x)

        return fe1
    else:
        raise NotImplementedError('bad input')


def Fanno_propagate(fe: FlowElement, friction, distance) -> FlowElement:

    xfinal = distance + fe.x

    M1 = Fanno_Mach_integrator(fe.M, xfinal, fe.Dh, fe.gamma, friction, x_initial=fe.x)
    p1 = Fanno_pressure_ratio(M1, fe.M, fe.gamma) * fe.p
    T1 = Fanno_temperature_ratio(M1, fe.M, fe.gamma) * fe.T

    fe1 = FlowElement(M1, p1, T1, fe.gamma, fe.R, fe.A, fe.x)

    return fe1

def shock(fe:FlowElement) -> FlowElement:
    pass
        
    





