import numpy as np
import math
from src.helper import *
from src.FlowElement import *

gamma = 1.4

friction = 0.0025

p0tot = 10 * 101325  # Pa
T0tot = 300

A0 = 1.0 # reference area

Atest__Astar = 3.0
L__D_test = 12.0

FE0 = FlowElement(0.0, p0tot, T0tot, A=A0)
FEstar = isentropic_propagate(FE0, M1=1.0)

Astar = FEstar.A
Atest = Atest__Astar * Astar
Dtest = math.sqrt(Atest * 4 / math.pi)
Ltest = Dtest * L__D_test

# a)
#1) shock at the throat
FEstar_prime1 = normal_shock(FEstar)
FEtest_section1 = isentropic_propagate(FEstar_prime1, A1=Atest)
FEexhaust1 = Fanno_propagate(FEtest_section1, friction=friction, distance=Ltest)
# print(FE0, FEstar, FEstar_prime1, FEtest_section1, FEexhaust1)
print('pressure ratio 1: {:.2f}'.format(FE0.p / FEexhaust1.p))
print('total pressure ratio 1: {:.2f}'.format(FE0.ptot / FEexhaust1.ptot))
print('exit pressure 1: {:.2f}'.format(FEexhaust1.p))

# 2 shock at nozzle exit
FEtest_section2 = isentropic_propagate(FEstar, A1=Atest)
FEtest_section2_shock = normal_shock(FEtest_section2)
FEexhaust2 = Fanno_propagate(FEtest_section2_shock, friction=friction, distance=Ltest)
# print(FE0, FEstar, FEtest_section2, FEtest_section2_shock, FEexhaust2)
print('pressure ratio 2: {:.2f}'.format(FE0.p / FEexhaust2.p))
print('total pressure ratio 2: {:.2f}'.format(FE0.ptot / FEexhaust2.ptot))
print('exit pressure 2: {:.2f}'.format(FEexhaust2.p))

#3) shock at duct exit
FEtest_section3 = isentropic_propagate(FEstar, A1=Atest)
FEexhaust3 = Fanno_propagate(FEtest_section3, friction=friction, distance=Ltest)
FEexhaust3_shock = normal_shock(FEexhaust3)
print(FE0, FEstar, FEtest_section2, FEexhaust3, FEexhaust3_shock)
print('pressure ratio 3: {:.2f}'.format(FE0.p / FEexhaust3_shock.p))
print('total pressure ratio 3: {:.2f}'.format(FE0.ptot / FEexhaust3_shock.ptot))
print('exit pressure 3: {:.2f}'.format(FEexhaust3_shock.p))

# b) no shocks!
FEtest_section_b = isentropic_propagate(FEstar, A1=Atest)
FEexhaust_b = Fanno_propagate(FEtest_section3, friction=friction, distance=Ltest)
print('pressure ratio b): {:.2f}'.format(FE0.p / FEexhaust_b.p))
print('total pressure ratio b): {:.2f}'.format(FE0.ptot / FEexhaust_b.ptot))
print('exit pressure b): {:.2f}'.format(FEexhaust_b.p))

# print(FE0, FEstar, FEtest_section_b, FEexhaust_b)


