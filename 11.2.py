import numpy as np
import math
from src.FlowElement import *
from scipy.optimize import minimize

gamma = 1.4

friction = 0.005
Dh = 0.05 # m
A = Dh**2 * math.pi / 4
L = 0.6 # m

p0 = 80 * 1000  # Pa
T0 = 20 + 273.15  # K
M0 = 2.0

FE0 = FlowElement(M0, p0, T0, A=A)
FE1 = Fanno_propagate(FE0, friction=friction, distance=L)
FEstar = isentropic_propagate(FE1, M1=1)
A_exit = FEstar.A * 3
FE2_from_1 = isentropic_propagate(FE1, A1=A_exit)
FE2_from_star = isentropic_propagate(FEstar, A1=A_exit)

print(FE1)
print(FEstar)
print(FE2_from_1)
print(FE2_from_star)

# p0__psonic = Fanno_p_over_p_sonic(M0, gamma)
# T0__Tsonic = Fanno_T_over_T_sonic(M0, gamma)
#
# M1 = Fanno_Mach_integrator(M0=M0, Dh=Dh, x_final=L, gamma=gamma, friction=friction)
#
# p1__psonic = Fanno_p_over_p_sonic(M1, gamma)
# T1__Tsonic = Fanno_T_over_T_sonic(M1, gamma)
#
# p1 = p1__psonic / p0__psonic * p0
# T1 = T1__Tsonic / T0__Tsonic * T0
#
# print("conditions at the pipe exit:")
# print("p1: {:.2f} Pa".format(p1))
# print("T1: {:.2f} K".format(T1))
# print("M1: {:.2f}".format(M1))
#
# # 2 - exhaust, star - sonic flow
#
# p1__ptot = isentropic_p_over_p_tot(M1, gamma)
# T1__tot = isentropic_T_over_T_tot(M1, gamma)
# A1__Astar = isentropic_A_over_Astar(M1, gamma)
# A1 = Dh**2 / 4 * math.pi # area of the nozzle inlet (pipe)
#
# # Astar = A1 / A1__Astar
# # A2 = 3 * Astar
# A2__Astar = 3.0
#
# M2 = isentropic_Mach_from_area_ratio(A2__Astar, gamma)
#
# p2__ptot = isentropic_p_over_p_tot(M2, gamma)
# p2 = p2__ptot / p1__ptot * p1
#
# print("conditions at the nozzle exit:")
# print("p2: {:.2f} Pa".format(p2))
# print("M2: {:.2f}".format(M2))
