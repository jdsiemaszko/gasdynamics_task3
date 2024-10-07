import numpy as np
import math
from src.FlowElement import *
gamma = 1.4

friction = 0.005
Dh = 0.05
A = Dh**2 * math.pi / 4

p0 = 250 * 1000  # Pa
T0 = 30 + 273.15  # K
M0 = 0.0

M1 = 0.3
M2 = 0.95

FE0 = FlowElement(M0, p0, T0, A=A)
FE1 = isentropic_propagate(FE0, M1=M1)
FE2 = Fanno_propagate(FE1, friction=friction, M1=M2)

L_pipe = FE2.x - FE1.x

L_ref = L_pipe * 0.75
FE3 = Fanno_propagate(FE1, friction=friction, distance=L_ref)

for fp in [FE0, FE1, FE2, FE3]:
    print(fp)


# print(L_pipe)

# def isentropic_expr_pressure(M, gamma = 1.4):
#     return (1 + (gamma-1)/2 * M**2) * (gamma / (gamma-1))
#
# p1 = p0 * isentropic_expr_pressure(M0, gamma) / isentropic_expr_pressure(M1, gamma)
#
# p1__psonic = Fanno_p_over_p_sonic(M1, gamma)
# p2__psonic = Fanno_p_over_p_sonic(M2, gamma)
# p2__p1 = p2__psonic / p1__psonic
# p2 = p1 * p2__p1
#
# L_pipe = Fanno_length_from_Mach_difference(M0=M1, Mfinal=M2, gamma=gamma, friction=friction, Dh=Dh)
#
# print("Pipe Length: {:.2f} m".format(L_pipe))
# print("Pressure: {:.2f} Pa".format(p2))
#
# L_ref = L_pipe * 0.75
# M_ref = Fanno_Mach_integrator(M0=M1, Dh=Dh, gamma=gamma, x_final=L_ref)
# p_ref__psonic = Fanno_p_over_p_sonic(M_ref, gamma)
# p_ref__p1 = p_ref__psonic / p1__psonic
# p_ref = p_ref__p1 * p1
#
# print("Mach at 75% length: {:.2f}".format(M_ref))
# print("Pressure: {:.2f} Pa".format(p_ref))

