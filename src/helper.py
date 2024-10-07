import numpy as np
from scipy.optimize import minimize
import math

def Fanno_Mach_integrator(M0, x_final, Dh, gamma=1.4, friction=3e-3, N=1000, x_initial=0.0):
    M_current = M0
    dx = (x_final - x_initial) / N

    for i in range(N):
        M_current += dx * (
                -2 * M_current ** 2 * (1 + (gamma - 1) / 2 * M_current ** 2) * gamma * friction /
                Dh * (M_current ** 2 - 1)
        )

    return M_current


def Fanno_length_from_Mach_difference(M0, Mfinal, Dh, x_initial=0.0, gamma=1.4, friction=3e-3, N=1000):
    res = minimize(lambda x: abs(Fanno_Mach_integrator(M0 = M0, x_final=x, Dh=Dh,
                                                       gamma=gamma, friction=friction, N=N, x_initial=x_initial) - Mfinal),
                   x0=Dh * 10,
                   bounds=[(0.0, 1000*Dh)],
                   method='Nelder-Mead'
                   )
    return res.x[0]

def Fanno_p_over_p_sonic(M, gamma=1.4):

    sqr = math.sqrt(
        (2 / (gamma+1)) * (1 + (gamma-1) / 2 * M**2)
    )
    return 1 / M / sqr

def Fanno_pressure_ratio(M0, M1, gamma=1.4):
    return Fanno_p_over_p_sonic(M0, gamma) / Fanno_p_over_p_sonic(M1, gamma)

def Fanno_T_over_T_sonic(M, gamma=1.4):
    return 1 / ((2 / (gamma+1)) * (1 + (gamma-1) / 2 * M**2))

def Fanno_temperature_ratio(M0, M1, gamma=1.4):
    return Fanno_T_over_T_sonic(M0, gamma) / Fanno_T_over_T_sonic(M1, gamma)

def isentropic_p_over_p_tot(M, gamma=1.4):
    return 1 / (1 + (gamma - 1) / 2 * M ** 2) * (gamma / (gamma - 1))

def isentropic_T_over_T_tot(M, gamma=1.4):
    return 1 / (1 + (gamma - 1) / 2 * M ** 2)

def isentropic_A_over_Astar(M, gamma=1.4):

    a = ((gamma+1) / 2) ** (-(gamma+1) / 2 / (gamma-1))
    b = (1+ (gamma-1)/2 * M**2) ** ((gamma+1) / 2 / (gamma-1))

    return a * b / M

def isentropic_Mach_from_area_ratio(A_Astar, gamma=1.4):
    # invert the isentropic_A_over_Astar function
    res = minimize(lambda Mach: abs(isentropic_A_over_Astar(Mach, gamma) - A_Astar),
                   x0=0.9,
                   bounds=[(0.01, 10.0)],
                   method='Nelder-Mead'
                   )

    M = res.x[0]
    return M

if __name__ == "__main__":
    Dh = 0.05**2
    M0 = 0.3
    for x_final in np.linspace(0.1, 1.0, 5):
        # print(x_final)
        m = Fanno_Mach_integrator(M0, x_final, Dh)
        # print(m)
        xpred = Fanno_length_from_Mach_difference(M0, m, Dh)
        # print(xpred)

        print(np.abs(x_final - xpred) / x_final * 100)




