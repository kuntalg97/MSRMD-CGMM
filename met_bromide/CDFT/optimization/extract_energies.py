from fit import *
from sys import argv

eg, h11, h22, h12 = read_evb(argv[1])

np.savetxt("energy.dat", np.array([eg, h11, h22, h12]).transpose(), header = 'FIELDS Eg E1 E2 V12', comments = '#! ')
