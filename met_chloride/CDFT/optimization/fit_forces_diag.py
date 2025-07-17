# Kuntal Ghosh (2023-24)

# Based on the code written by Chenghan Li (J. Phys. Chem. B 2021, 125, 37, 10471–10480)

from evb_out import *
from read_colvar import *
import re
from os import system, path
from scipy.optimize import minimize
from scipy.linalg import sqrtm

# off diag to be used
offdiag_name = "METCL_SN2_PT"

# ref DFT energies
df_ref_name = "ci_ts.dat"
hartree = 627.5096

# initial parameters of Nelson 
# Lorentz-Berthelot rules for initial mixing parameters
params = dict()
params['global_morse_cut']    = 3.0
params['De_ccl_intra']        = 2.0
params['De_ccl_inter']        = 2.0
params['alpha_ccl_intra']     = 1.5
params['alpha_ccl_inter']     = 1.5
params['re_ccl_morse']        = 4.3
params['OFFDIAG_G']           = 5.0 
params['V_IJ_CONSTANT']       = -50.0
params['DIAG_VII_CONST']      = 0.000
names = ['global_morse_cut', 'De_ccl_intra', 'De_ccl_inter', 'alpha_ccl_intra', 'alpha_ccl_inter', 're_ccl_morse', 'OFFDIAG_G', 'V_IJ_CONSTANT', 'DIAG_VII_CONST']

bounds = np.array([[0,None],[None,None],[None,None],[0,None],[0,None],[0,None],[0,None],[None,None],[None,None]])

def write_input(params):
    fpin = open("sn2_template.inp")
    fpout = open("sn2_interim.inp", "w")

    for line in fpin:
       line = re.sub("global_morse_cut", str(params['global_morse_cut']), line)
       line = re.sub("De_ccl_intra", str(params['De_ccl_intra']), line)
       line = re.sub("De_ccl_inter", str(params['De_ccl_inter']), line)
       line = re.sub("alpha_ccl_intra", str(params['alpha_ccl_intra']), line)
       line = re.sub("alpha_ccl_inter", str(params['alpha_ccl_inter']), line)
       line = re.sub("re_ccl_morse", str(params['re_ccl_morse']), line)
       fpout.write(line)

    fpin.close()
    fpout.close()

    fpout = open("evb_model.par", "w")
    print("#define", offdiag_name, file = fpout)
    for i in range(len(names)):
       n = names[i]
       if n in ['global_morse_cut', 'De_ccl_intra', 'De_ccl_inter', 'alpha_ccl_intra', 'alpha_ccl_inter', 're_ccl_morse']:
          continue
       print("#define {:15s} {:<20.15f}".format(n, params[n]), file = fpout)

    fpout.close()

def array2params(arr):
    params = dict()
    for i in range(len(names)):
       params[names[i]] =  arr[i]
    return params

def params2array(params):
    arr = list()
    for i in range(len(names)):
       arr.append(params[names[i]])
    return arr

def get_forces(dist, energies):
    npoints = np.size(energies)
    dx = dist[1]-dist[0]

    forces = np.empty(npoints)
    for i in range (npoints-1):
        forces[i] = -(energies[i+1]-energies[i])/dx

    forces = forces [1:-1]
    return forces

def loss(params):
    write_input(params)
    if path.exists("evb.out"):
       system("rm evb.out")
    system("lmp -in sn2_interim.inp > out")
#    print ("New LAMMPS run")
    Eg, E1, E2, V12 = read_evb("evb.out")
    Eg_ref = df_ref['Eg'].values * hartree 
    V12_ref = df_ref['V12'].values * hartree
    E1_ref = df_ref['E1'].values * hartree
    E2_ref = df_ref['E2'].values * hartree
    c1_ref = df_ref['normed_c1']
    c2_ref = df_ref['normed_c2']
    dist = df_ref['dist_diff']
    offset = min(Eg) - min(Eg_ref)

    forces_Eg_ref = get_forces(dist, Eg_ref)
    forces_Eg = get_forces(dist, Eg)

    Hs = np.zeros((len(Eg),2,2))
    Hs[:,0,0] = E1
    Hs[:,1,1] = E2
    Hs[:,0,1] = V12
    Hs[:,1,0] = V12

#   Orthogonalizing the Hamiltonian

    S12_ref = df_ref['S12']
    S12 = np.zeros((len(Eg),2,2))
    Hs_ortho = np.zeros((len(Eg),2,2))
    S12[:,0,0] = 1
    S12[:,1,1] = 1
    S12[:,0,1] = S12_ref
    S12[:,1,0] = S12_ref

    for i in  range(len(Eg)):
        S12_sqrt_inv = np.linalg.inv(sqrtm(S12[i]))
        Hs_ortho[i] = np.dot(np.dot(S12_sqrt_inv, Hs[i]), S12_sqrt_inv)

#    _, vectors = np.linalg.eigh(Hs)    # eigenvectors already normalized
    _, vectors = np.linalg.eigh(Hs_ortho)    # eigenvectors already normalized
    ee  = np.mean((abs(vectors[:,0,0]) - c1_ref)**2) 
    ee += np.mean((abs(vectors[:,1,0]) - c2_ref)**2)
    return ee

def loss_arr(arr):
    global call_counter
    print ("iter", call_counter)
    global fplog
    call_counter += 1
    if call_counter % 10 == 0:
        p = array2params(arr)
        print(p, loss(p), file = fplog, flush = True)
    return loss(array2params(arr)) 

def save_params(arr, fname):
    params = array2params(arr)
    fp = open(fname, "w")
    for i in range(len(names)):
       n = names[i]
       print("{:<20.15f} # {:s}".format(params[n],n), file = fp)
    fp.close()

def read_evb(fname, get_ci = False):
    evb = evb_out(fname)
    env = evb.ene_environment_trj
    eg = evb.ene_total_trj
    h11 = []; h22 = []; h12 = []
    if get_ci:
        c1 = []
        c2 = []
    for i in range(len(env)):
        diag = evb.diag_trj[i]
        states = evb.states_trj[i]
        offdiag = evb.offdiag_trj[i]
        eigenv = evb.eigenv_trj[i]

        sel_glup = np.where(states.mol_B == 1)[0]
        glup_istate = sel_glup[0]
        glup_ene = diag.total[glup_istate]

        sel_h3o = np.where(states.mol_B == 2)[0]
        h3o_istate = sel_h3o[0]
        h3o_ene = diag.total[h3o_istate]

        offdiag_ene = offdiag.energy[0]

        h11.append(glup_ene)
        h22.append(h3o_ene)
        h12.append(offdiag_ene)

        if get_ci:
            c1.append(eigenv[glup_istate])
            c2.append(eigenv[h3o_istate])

    if get_ci:
       return [np.array(eg), np.array(h11)+env, np.array(h22)+env, np.array(h12), \
               np.array(c1), np.array(c2)]
    else:
       return [np.array(eg), np.array(h11)+env, np.array(h22)+env, np.array(h12)]

if __name__ == "__main__":
    global call_counter
    call_counter = 1
    global fplog
    global df_ref
    df_ref = read_colvar(df_ref_name)
    method = "Nelder-Mead"
    #method = "BFGS"
    fplog = open(f"fit_{method}.log", "w")
    init = params2array(params)
    loss_arr(init)  
    opt = minimize(loss_arr, init, method = method, options = {'fatol': 1e-5}, bounds = bounds)
    #opt = minimize(loss_arr, init, method = method, options = {'fatol': 1e-5})
    #opt = minimize(loss_arr, init, method = method, options = {'gtol': 1e-5})
    save_params(opt['x'], f"opt_{method}.dat")
    fplog.close()
