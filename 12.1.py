import numpy as np
import math
from src.helper import *
from src.FlowElement import *

FAR = 1/40 # fuel to air ratio
p0 = 50 * 1000
T0 = 30 + 273.15
V0 = 80
gamma = 1.4
R = 287.05
Cp = R * gamma / (gamma-1)

A_ref = 1.0

# heating value
alpha = 40 * 10**6 # J/kg

a0 = math.sqrt(gamma * R * T0)
M0 = V0 / a0
FP0 = FlowElement(M0, p0, T0, A=A_ref)

delta_q = alpha * FAR * V0 * A_ref * FP0.rho # total amount of heat added per unit time

FP1 = Rayleigh_heat_addition(FP0, delta_q) # propagate along Rayleigh line

print(FP0, FP1)


