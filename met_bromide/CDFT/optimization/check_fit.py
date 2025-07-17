#from colors import *
import matplotlib.pyplot as plt
import numpy as np

from sys import argv
from sys import path
from read_colvar import *
c_ref = read_colvar("ci_ts.dat")
c_evb = read_colvar(argv[1])

fig, ax = plt.subplots()
ax.set_xlim((0,1))
ax.set_ylim((0,1))
ax.set_xlabel(r"$c^\mathrm{RMD}$")
ax.set_ylabel(r"$c^\mathrm{CDFT}_\mathrm{normed}$")
ax.scatter(list(c_evb['c1'])+list(c_evb['c2']), list(c_ref['normed_c1'])+list(c_ref['normed_c2']), s=50)
#ax.scatter(list(c_evb['c1']**2)+list(c_evb['c2']**2), list(c_ref['normed_c1']**2)+list(c_ref['normed_c2']**2), s=200)
#ax.scatter(list(c_evb['c2']), list(c_ref['normed_c2']), marker='x', s=200)
#ax.plot((0,1),(0,1),c=tableau_10[1])
ax.plot((0,1),(0,1))
plt.savefig(f"scatterplot_{argv[1][:-4]}.png", bbox_inches='tight')

print("<|c^EVB-c^CDFT|^2> =", sum(np.mean((c_evb[['c1','c2']].values - c_ref[['normed_c1','normed_c2']])**2)))
