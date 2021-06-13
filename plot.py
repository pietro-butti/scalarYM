import numpy as np
import math
import os
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from operator import itemgetter


ciao = np.loadtxt('prova.dat')


plt.plot(ciao)
plt.show()