from fit import *
from sys import argv

eg, h11, h22, h12, c1, c2 = read_evb(argv[1], get_ci = True)

np.savetxt("ci.dat", np.array([c1, c2]).transpose(), header = 'FIELDS c1 c2', comments = '#! ')
