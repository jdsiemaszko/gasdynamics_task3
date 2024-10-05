import numpy as np
import math
from src.helper import *

gamma = 1.4

friction = 0.005
Dh = 0.05 # m
L = 0.6 # m

p0 = 80 * 1000  # Pa
T0 = 20 + 273.15  # K
M0 = 2.0

p0__psonic = Fanno_p_over_p_sonic(M0, gamma)
T0__Tsonic = Fanno_T_over_T_sonic(M0, gamma)


M1 = Fanno_Mach_integrator(M0=M0, Dh=Dh, x_final=L, gamma=gamma, friction=friction)

p1__psonic = Fanno_p_over_p_sonic(M1, gamma)
T1__Tsonic = Fanno_T_over_T_sonic(M1, gamma)

p1 = p1__psonic / p0__psonic * p0
T1 = T1__Tsonic / T0__Tsonic * T0

print("conditions at the pipe exit:")
print("p1: {:.2f} Pa".format(p1))
print("T1: {:.2f} K".format(T1))
print("M1: {:.2f}".format(M1))

# 2 - exhaust, star - sonic flow

p1__ptot = isentropic_p_over_p_tot(M1, gamma)
T1__tot = isentropic_T_over_T_tot(M1, gamma)


