import numpy as np
import math
from src.helper import *

gamma = 1.4

friction = 0.0025
Dh = 0.05 # m
L = 0.6 # m

p0tot = 10 * 101325  # Pa
T0tot = 300

A1__Astar = 3.0

Dh = 1.0
L = 12 * Dh

#1) shock at the throat
# assume very high pressure difference!

