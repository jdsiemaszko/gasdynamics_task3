"""
Helper functions for Task 3, including Fanno, Rayleigh, and isentropic relations
"""

import numpy as np
from scipy.optimize import minimize
import math

def Fanno_integrator(M0, p0, T0, x_final, Dh, gamma=1.4, friction=3e-3, N=1000, x_initial=0.0):
    """
    manually inegrate the Mach number along the Fanno line
    note: p and t can be computed as functions of the Mach number, so the results computed here are not used
    """

    M_current = M0
    T_current = T0
    p_current = p0

    dx = (x_final - x_initial) / N

    for i in range(N):
        M_current += dx * (
                -2 * M_current ** 2 * (1 + (gamma - 1) / 2 * M_current ** 2) * gamma * friction /
                Dh * (M_current ** 2 - 1)
        )
        if M_current != 1:
            p_current += dx * friction * Dh * ((-gamma*M_current**2*(1+(gamma-1)*M_current**2))/(1-M_current**2)) * p_current
            T_current += dx * friction * Dh * ((-gamma * (gamma-1)*M_current**4)/(1-M_current**2)) * T_current

    return M_current, p_current, T_current


def Fanno_length_from_Mach_difference(M0, Mfinal, Dh, x_initial=0.0, gamma=1.4, friction=3e-3, N=1000):
    """
    inverse of Fanno_integrator to find pipe length from Minitial and Mfinal
    """
    res = minimize(lambda x: abs(Fanno_integrator(M0, 0.0, 0.0, x_final=x, Dh=Dh,
                                                       gamma=gamma, friction=friction, N=N, x_initial=x_initial) - Mfinal),
                   x0=Dh * 10,
                   bounds=[(0.0, 1000*Dh)],
                   method='Nelder-Mead'
                   )
    return res.x[0]

def Fanno_max_unchoked_length(M0, gamma=1.4, Dhref = 1.0, fref = 1.0):
    # as per the reader
    a = (1 - M0**2) / M0**2
    b = (gamma+1)/2 * math.log(
        (gamma+1) * M0**2 / (2*(1+(gamma-1)/2*M0**2))
    )

    Lmax = (a+b) * Dhref / 4 / gamma / fref
    return Lmax

def Fanno_M0_from_choked_length(Lchoked, gamma=1.4, Dhref = 1.0, fref = 1.0):
    # inverse of Fanno_max_unchoked_length to find M0 from Lmax
    res = minimize(lambda M0: abs(Fanno_max_unchoked_length(M0, gamma, Dhref, fref) - Lchoked),
                   x0=0.2,
                   bounds=[(0.0, 1.0)], # always subsonic?
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
    return 1 / (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))

def isentropic_pressure_ratio(M0, M1, gamma=1.4):
    return isentropic_p_over_p_tot(M0, gamma) / isentropic_p_over_p_tot(M1, gamma)

def isentropic_T_over_T_tot(M, gamma=1.4):
    return 1 / (1 + (gamma - 1) / 2 * M ** 2)

def isentropic_temp_ratio(M0, M1, gamma=1.4):
    return isentropic_T_over_T_tot(M0, gamma) / isentropic_T_over_T_tot(M1, gamma)

def isentropic_A_over_Astar(M, gamma=1.4):
    if M == 0:
        return 0.0

    a = ((gamma+1) / 2) ** (-(gamma+1) / 2 / (gamma-1))
    b = (1+ (gamma-1)/2 * M**2) ** ((gamma+1) / 2 / (gamma-1))

    return a * b / M

def isentropic_area_ratio(M0, M1, gamma=1.4):
    if any([M1 == 0, M0 == 0]):
        return None
    return isentropic_A_over_Astar(M0, gamma) / isentropic_A_over_Astar(M1, gamma)

def isentropic_Mach_from_area_ratio(A1_A0, M0, gamma=1.4):

    # bounds for the inverse function
    a = 1e-6 if M0 < 1 else 1 + 1e-6
    b = 1 - 1e-6 if M0 < 1 else 20
    x0 = 0.5 if M0 < 1 else 1.5

    res = minimize(lambda Mach: abs(isentropic_area_ratio(Mach, M0, gamma) - A1_A0),
                   x0=x0,
                   bounds=[(a, b)],
                   method='Nelder-Mead'
                   )

    M1 = res.x[0]
    return M1

def normal_shock_upstream_mach(M2, gamma=1.4):
    res = minimize(lambda M1: abs((math.sqrt(((gamma - 1) * M1**2 + 2) / (2 * gamma * M1**2 - (gamma - 1))) if M1 > 1 else 1-1e-6) - M2),
                   x0=1.5,
                   bounds=[(1+1e-6, 10)],
                   method='Nelder-Mead'
                   )

    M1 = res.x[0]
    return M1

def Rayleigh_T_over_Tstar(M, gamma=1.4):
    return (1+gamma**2)*M**2 / (1+gamma*M**2)**2

def Rayleigh_T_ratio(M1, M2, gamma=1.4):
    return Rayleigh_T_over_Tstar(M1, gamma) /  Rayleigh_T_over_Tstar(M2, gamma)

def Rayleigh_p_over_pstar(M, gamma=1.4):
    return (1+gamma) / (1+gamma*M**2)**2

def Rayleigh_p_ratio(M1, M2, gamma=1.4):
    return Rayleigh_p_over_pstar(M1, gamma) / Rayleigh_p_over_pstar(M1, gamma)

def __Rayleigh_Ttot_helper(M, gamma=1.4):
    return M**2 * (1+(gamma-1)*M**2/2)/(1+gamma*M**2)

def Rayleigh_Ttot_ratio_from_Mach(M1, M2, gamma=1.4):
    a = __Rayleigh_Ttot_helper(M1, gamma)
    b = __Rayleigh_Ttot_helper(M2, gamma)
    return a/b

def Rayleigh_Mach_from_Ttot_ratio(Ttot1, Ttot2, M1, gamma=1.4):

    a = 1e-6 if M1 < 1 else 1 + 1e-6
    b = 1 - 1e-6 if M1 < 1 else 20
    x0 = 0.5 if M1 < 1 else 1.5

    val = Ttot1 / Ttot2
    res = minimize(lambda M2: abs(Rayleigh_Ttot_ratio_from_Mach(M2, M1, gamma) - val),
                   x0=x0,
                   bounds=[(a, b)],
                   method='Nelder-Mead'
                   )

    M2 = res.x[0]
    return M2









