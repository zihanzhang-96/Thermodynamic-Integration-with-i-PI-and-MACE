import numpy as np
import os, glob
from ase.io import read, write
from collections import Counter
from ase.data import atomic_masses, atomic_numbers

#data_pattern = "*.vasp"
#p_path = '/home/users/nus/zhongpc/scratch/zzh/lmp-mace/MACE-matpes-r2scan-omat-ft.model-lammps.pt'
#files = sorted(glob.glob(data_pattern))
#T = 100
#L = 0.3

def input(pos, p_path, L, T):
    ff = pos.split('.')[0]
    L_T_dir = f"{L:.4f}"+'_'+f"{T:.1f}"+'_dir'
    ff2 = ff+'/I_L/'+L_T_dir

    atoms = read(ff+"/POSCAR")
    symbols = atoms.get_chemical_symbols()
    count = Counter(symbols)

    elem = []
    num = []
    mass = []
    for e in count:
        elem.append(e)
        num.append(count[e])
        Z = atomic_numbers[e]
        mass.append(atomic_masses[Z])
    N = sum(num)

    os.system('mkdir '+ff+'/I_L ')
    os.system('mkdir '+ff2+' && cp '+ff+'/SUPERCELL.xyz '+ff+'/SUPERCELL.lmp '+ff+'/hes/perfect.perfect-fd.hess '+ff2)

    file_hes = open(ff2+'/input.xml','w')
    file_hes.write('<simulation verbosity="medium">'+'\n')
    file_hes.write("    <output prefix='perfect'>"+'\n')
    file_hes.write("        <properties stride='10' filename='out'>  [ step, time{picosecond}, conserved, temperature{kelvin}, potential{electronvolt}, pot_component_raw(0){electronvolt}, pot_component_raw(1){electronvolt}, ensemble_bias, ensemble_temperature{kelvin} ] </properties>"+'\n')
    file_hes.write("        <trajectory filename='pos' stride='1000' bead='0' format='xyz' cell_units='angstrom'> positions{angstrom} </trajectory>"+'\n')
    file_hes.write("    </output>"+'\n')
    file_hes.write("    <total_steps>200000</total_steps>"+'\n')
    file_hes.write("    <prng><seed>85556</seed></prng>"+'\n')
    file_hes.write('    <ffsocket name="lmpserial1" mode="unix" pbc="true">'+'\n')
    file_hes.write("        <address>localhost</address>"+'\n')
    file_hes.write("    </ffsocket>"+'\n')
    file_hes.write("    <ffdebye name='debye'>"+'\n')
    file_hes.write("        <hessian mode='file' shape='("+str(int(3*N))+","+str(int(3*N))+")'> perfect.perfect-fd.hess </hessian>"+'\n')
    file_hes.write("        <x_reference mode='file' units='angstrom'> x_reference.txt </x_reference>"+'\n')
    file_hes.write("    </ffdebye>"+'\n')
    file_hes.write("    <system>"+'\n')
    file_hes.write("        <initialize nbeads='1'>"+'\n')
    file_hes.write("            <file mode='xyz'> SUPERCELL.xyz </file>"+'\n')
    file_hes.write('            <velocities mode="thermal" units="kelvin"> '+str(T)+' </velocities>'+'\n')
    file_hes.write("        </initialize>"+'\n')
    file_hes.write("        <forces>"+'\n')
    file_hes.write('            <force forcefield="debye" weight="'+f"{1-L:.5f}"+'"> </force>'+'\n')
    file_hes.write('            <force forcefield="lmpserial1" weight="'+f"{L:.5f}"+'"> </force>'+'\n')
    file_hes.write("        </forces>"+'\n')
    file_hes.write('        <motion mode="dynamics">'+'\n')
    file_hes.write('            <dynamics mode="nvt">'+'\n')
    file_hes.write("                <thermostat mode='langevin'>"+'\n')
    file_hes.write("                    <tau units='femtosecond'>100</tau>"+'\n')
    file_hes.write("                </thermostat>"+'\n')
    file_hes.write("                <timestep units='femtosecond'> 1 </timestep>"+'\n')
    file_hes.write("            </dynamics> <fixcom> True </fixcom>"+'\n')
    file_hes.write("        </motion>"+'\n')
    file_hes.write("        <ensemble>"+'\n')
    file_hes.write("            <temperature units='kelvin'> "+str(T)+" </temperature>"+'\n')
    file_hes.write("            <pressure units='bar'> 1 </pressure>"+'\n')
    file_hes.write("        </ensemble>"+'\n')
    file_hes.write("    </system>"+'\n')
    file_hes.write("</simulation>"+'\n')
    file_hes.close()

    file_hes = open(ff2+'/in.lammps','w')
    file_hes.write('# I_L'+'\n'+'\n')
    file_hes.write("units           metal"+'\n')
    file_hes.write("boundary        p p p"+'\n')
    file_hes.write("atom_style      atomic"+'\n')
    file_hes.write("atom_modify     map yes"+'\n')
    file_hes.write("newton          on"+'\n'+'\n')
    file_hes.write("neighbor        2.0 bin"+'\n')
    file_hes.write("neigh_modify    every 10 delay 0 check no"+'\n'+'\n')
    file_hes.write("read_data       SUPERCELL.lmp"+'\n')
    for i in range(len(elem)):
        file_hes.write("mass            "+str(i+1)+f" {mass[i]:.3f}"+'\n')
    file_hes.write('\n'+'\n')
    file_hes.write("pair_style      mace no_domain_decomposition"+'\n')
    file_hes.write("pair_coeff      * * "+p_path+' ')
    for i in range(len(elem)):
        file_hes.write(elem[i]+" ")
    file_hes.write('\n'+'\n')
    file_hes.write("velocity        all create "+str(T)+" 8257 dist gaussian mom yes"+'\n'+'\n')
    file_hes.write("fix             1 all ipi localhost 31415 unix"+'\n')
    file_hes.write("timestep        0.0005"+'\n')
    file_hes.write("thermo_style    custom step pe ke etotal temp press vol"+'\n')
    file_hes.write("thermo          100"+'\n')
    file_hes.write("dump            1 all custom 10 SUPERCELL.dump id type x y z"+'\n'+'\n')
    file_hes.write("run             200000"+'\n')
    file_hes.close()

    file_hes = open(ff2+'/SUPERCELL.xyz','r')
    s = file_hes.readlines()
    file_hes.close()

    file_hes = open(ff2+'/x_reference.txt','w')
    for i in range(len(s)-2):
        file_hes.write(s[2+i][2:])
    file_hes.close()
