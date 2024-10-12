import numpy as np
import math
from src.helper import *
from src.FlowElement import *

T0 = 300
M0 = 1.5
p0 = 1e5 #reference

gamma = 1.4
M1 = 1.0
A_Astar = 1 / 0.98

FP0 = FlowElement(M0, p0, T0, gamma=gamma, R=287.05, A = 1.0)

M_inlet = isentropic_Mach_from_area_ratio(A_Astar, 1.0-1e-6, gamma)
M_preshock = normal_shock_upstream_mach(M_inlet, gamma) # invert shock relations to find upstream Mach

FP1 = Rayleigh_propagate(FP0, M_preshock)
FP1_shock = normal_shock(FP1)
print(M_inlet, FP1_shock.M) # check if Mach matches what's expected
FP2 = isentropic_propagate(FP1_shock, A1 = FP1_shock.A / A_Astar)

q = (FP1.Ttot / FP0.Ttot - 1) * FP0.Cp * FP0.Ttot
print(q)

print(FP0, FP1, FP1_shock, FP2)



