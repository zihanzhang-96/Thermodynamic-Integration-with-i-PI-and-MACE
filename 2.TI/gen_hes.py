import numpy as np
import os, glob
from ase.io import read, write
from collections import Counter
from ase.data import atomic_masses, atomic_numbers

#data_pattern = "*.vasp"
#p_path = '/home/users/nus/zhongpc/scratch/zzh/lmp-mace/MACE-matpes-r2scan-omat-ft.model-lammps.pt'
#files = sorted(glob.glob(data_pattern))
#T = 100
#L = 1.0


def hes(pos, p_path):
    ff = pos.split('.')[0]

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

    os.system('mkdir '+ff+'/hes && cp '+ff+'/SUPERCELL.xyz '+ff+'/SUPERCELL.lmp '+ff+'/hes')   

    file_hes = open(ff+'/hes/input.xml','w')
    file_hes.write('<simulation verbosity="medium">'+'\n')
    file_hes.write("    <output prefix='perfect'>"+'\n')
    file_hes.write("        <properties stride='1' filename='out'>  [ step, time{picosecond}, conserved, temperature{kelvin}, potential, pressure_md, volume, cell_h ] </properties>"+'\n')
    file_hes.write("        <trajectory filename='pos' stride='1' bead='0' format='xyz' cell_units='angstrom'> positions{angstrom} </trajectory>"+'\n')
    file_hes.write("    </output>"+'\n')
    file_hes.write("    <total_steps>100000</total_steps>"+'\n')
    file_hes.write("    <prng><seed>832346</seed></prng>"+'\n')
    file_hes.write('    <ffsocket name="lmpserial1" mode="unix" pbc="true">'+'\n')
    file_hes.write("        <address>localhost</address>"+'\n')
    file_hes.write("    </ffsocket>"+'\n')
    file_hes.write("    <system>"+'\n')
    file_hes.write("        <initialize nbeads='1'>"+'\n')
    file_hes.write("            <file mode='xyz'> SUPERCELL.xyz </file>"+'\n')
    file_hes.write('            <velocities mode="thermal" units="kelvin"> 300 </velocities>'+'\n')
    file_hes.write("        </initialize>"+'\n')
    file_hes.write("        <forces>"+'\n')
    file_hes.write('            <force forcefield="lmpserial1"> </force>'+'\n')
    file_hes.write("        </forces>"+'\n')
    file_hes.write('        <motion mode="vibrations">'+'\n')
    file_hes.write('            <vibrations mode="fd">'+'\n')
    file_hes.write("                <pos_shift> 0.001 </pos_shift>"+'\n')
    file_hes.write("                <energy_shift> 0.0009500 </energy_shift>"+'\n')
    file_hes.write("                <prefix> perfect-fd </prefix>"+'\n')
    file_hes.write("                <asr> crystal </asr>"+'\n')
    file_hes.write("            </vibrations>"+'\n')
    file_hes.write("        </motion>"+'\n')
    file_hes.write("        <ensemble>"+'\n')
    file_hes.write("        <temperature units='kelvin'> 300 </temperature>"+'\n')
    file_hes.write("        <pressure units='megapascal'> 1e-10 </pressure>"+'\n')
    file_hes.write("        </ensemble>"+'\n')
    file_hes.write("    </system>"+'\n')
    file_hes.write("</simulation>"+'\n')
    file_hes.close()

    file_hes = open(ff+'/hes/in.lammps','w')
    file_hes.write('# hes'+'\n'+'\n')
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
    file_hes.write("velocity        all create 300.0 8257 dist gaussian mom yes"+'\n'+'\n')
    file_hes.write("fix             1 all ipi localhost 31415 unix"+'\n')
    file_hes.write("timestep        0.0005"+'\n')
    file_hes.write("thermo_style    custom step pe ke etotal temp press vol"+'\n')
    file_hes.write("thermo          100"+'\n')
    file_hes.write("dump            1 all custom 10 SUPERCELL.dump id type x y z"+'\n'+'\n')
    file_hes.write("run             10000"+'\n')
    file_hes.close()
