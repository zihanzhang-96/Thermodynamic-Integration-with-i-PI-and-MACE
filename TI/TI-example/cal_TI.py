import numpy as np
import os, glob

import gen_sc, gen_hes, gen_input

data_pattern = "*.vasp"
p_path = '/home/users/nus/zhongpc/scratch/zzh/lmp-mace/MACE-matpes-r2scan-omat-ft.model-lammps.pt'
run_lmp = 'mpirun /home/users/nus/zhongpc/softwares/mace-lmp-ipi/lammps/build-ampere/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.lammps'
L = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
T = [100]

L2 = [1.0]
T2 = [100,  117,  138,  162,  190,  224,  263,  309,  362,  425,  500]

files = sorted(glob.glob(data_pattern))

for pos in files:
    ff = pos.split('.')[0]
#####################################supercell############################
    gen_sc.SC(pos)
#####################################calculate hessian####################
    gen_hes.hes(pos, p_path)
    os.system('rm /tmp/ipi_localhost')
    os.system('cd '+ff+'/hes && i-pi input.xml &')
    os.system('sleep 20')
    if os.path.isfile(ff+"/hes/perfect.out"):
        os.system('cd '+ff+'/hes && '+run_lmp)
    else:
        os.system('sleep 10')
    if os.path.isfile(ff+"/hes/perfect.out"):
        os.system('cd '+ff+'/hes && '+run_lmp)
    else:
        os.system('sleep 10')
    os.system('cd '+ff+'/hes && '+run_lmp)
#####################################Integral with L######################
    for i in T:
        for j in L:
            L_T_dir = f"{j:.4f}"+'_'+f"{i:.1f}"+'_dir'
            ff2 = ff+'/I_L/'+L_T_dir
            gen_input.input(pos, p_path, j, i)
            os.system('rm /tmp/ipi_localhost')
            os.system('cd '+ff2+' && i-pi input.xml &')
            os.system('sleep 20')
            if os.path.isfile(ff2+"/perfect.out"):
                os.system('cd '+ff2+' && '+run_lmp)
            else:
                os.system('sleep 10')
            if os.path.isfile(ff2+"/perfect.out"):
                os.system('cd '+ff2+' && '+run_lmp)
            else:
                os.system('sleep 10')
            os.system('cd '+ff2+' && '+run_lmp)
#####################################Integral with T######################
    for i in T2:
        for j in L2:
            L_T_dir = f"{j:.4f}"+'_'+f"{i:.1f}"+'_dir'
            ff2 = ff+'/I_L/'+L_T_dir
            gen_input.input(pos, p_path, j, i)
            os.system('rm /tmp/ipi_localhost')
            os.system('cd '+ff2+' && i-pi input.xml &')
            os.system('sleep 20')
            if os.path.isfile(ff2+"/perfect.out"):
                os.system('cd '+ff2+' && '+run_lmp)
            else:
                os.system('sleep 10')
            if os.path.isfile(ff2+"/perfect.out"):
                os.system('cd '+ff2+' && '+run_lmp)
            else:
                os.system('sleep 10')
            os.system('cd '+ff2+' && '+run_lmp)
