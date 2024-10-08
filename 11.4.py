import numpy as np
import math
from src.helper import *
from src.FlowElement import *

gamma = 1.4

friction_length = 5.3

# reference values
friction_ref = 0.003
D_ref = 1.0
A_ref = D_ref**2 * math.pi / 4
L_ref = friction_length * D_ref / friction_ref

p0tot = 7 * 101325  # Pa
T0tot = 555 # K

FE1 = FlowElement(0.0, p0tot, T0tot, gamma=gamma, A = A_ref)
M2 = Fanno_M0_from_choked_length(L_ref, gamma, D_ref, friction_ref) # find the Mach at duct inlet
FE2 = isentropic_propagate(FE1, M1=M2)

FE4 = Fanno_propagate(FE2, friction=friction_ref, distance=L_ref)

p_rec = FE4.p
L_short = L_ref / 5

# shorter duct -> flow is now underexpanded at the exit!
FE2 = isentropic_propagate(FE1, M1=M2)  # point 2 unchanged
FE3_prime = Fanno_propagate(FE2, friction=friction_ref, distance=L_short)

# print(FE1, FE2, FE3, FE4)

print(FE3_prime.p, FE4.p)

