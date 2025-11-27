import numpy as np
import os, glob
import datetime
import random
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from ase.constraints import ExpCellFilter
from ase import units
#from ase.spacegroup import get_spacegroup
from collections import Counter
import spglib
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.file_IO import write_FORCE_SETS
from ase.vibrations import Vibrations

calc = MACECalculator(model_paths="/home/users/nus/zhongpc/zzh/mace_model/MACE-matpes-r2scan-omat-ft.model", device="cuda") 
delta = 0.01                       # finite difference displacement (Å)
temperatures = [0, 100, 200, 300, 400, 500]  # K 要计算的温度点
T_ref = 300

data_pattern = "*.vasp"

files = sorted(glob.glob(data_pattern))

os.system('mkdir G H opt phonon')

for pos in files:
    atoms = read(pos)
    atoms.set_pbc([True, True, True])
    atoms.set_calculator(calc)
    ecf = ExpCellFilter(atoms,mask=[1, 1, 1, 1, 1, 1])
    dyn = BFGS(ecf, trajectory=pos.split('.')[0] + ".traj", logfile=pos.split('.')[0] + ".log")
    dyn.run(fmax=0.01, steps=200)
    write("optimized_"+pos, atoms) 

    fileRes = open(pos.split('.')[0]+'.res','w')

    stress = atoms.get_stress(voigt=True)
    p = -np.mean(stress[:3]) * 160.21766208
    v = atoms.get_volume()
    e = atoms.get_potential_energy()
    TOT_Atoms = len(atoms)
    cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    try:
        dataset = spglib.get_symmetry_dataset(cell, symprec=1e-2)
        sg = dataset["international"]
    except Exception:
        sg = "Unknown"
    SYS = sg

    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    n_elements = len(counts)

    fileRes.write('TITL '+pos+'  '+str(p)+'  '+str(v)+'  '+str(e)+' 0 0 '+str(TOT_Atoms)+' ('+SYS+')  n - 1'+'\n')
    fileRes.write('REM'+'\n')
    #fileRes.write('REM Run started in '+ path +'\n')
    fileRes.write('REM'+'\n')
    fileRes.write('CELL 1.54180    '+str(round(a,6))+'    '+str(round(b,6))+'    '+str(round(c,6))+'    '+str(round(alpha,6))+'    '+str(round(beta,6))+'    '+str(round(gamma,6))+'    '+'\n')
    fileRes.write('LATT -1'+'\n')
    fileRes.write('SFAC ')
    for elem, num in counts.items():
        fileRes.write(elem+'  ')
    fileRes.write('\n')
    for j, atom in enumerate(atoms):
        x, y, z = atoms.get_scaled_positions()[j]
        fileRes.write(str(atom.symbol)+'  ')
        fileRes.write(str(1)+'  ')
        fileRes.write(str(x)+'  ')
        fileRes.write(str(y)+'  ')
        fileRes.write(str(z)+'  ')
        fileRes.write('1.0  ')
        fileRes.write('\n')

    fileRes.write('END'+'\n')
    fileRes.close()


    os.system('mv '+pos.split('.')[0]+'.res'+' H')
    os.system('cp '+"optimized_"+pos+' POSCAR')
    os.system('mv *traj *log '+"optimized_"+pos+' opt')


    file = open('vaspkit.in','w')
    file.write('305'+'\n')
    file.write('3'+'\n')
    file.close()
    os.system('vaspkit < vaspkit.in >>vaspkit.log')
    os.system('rm PRIMCELL.vasp')
    os.system('mkdir phonon/'+pos.split('.')[0])
    os.system("sed -i 's/^DIM = .*/DIM = 2 2 2/' KPATH.phonopy ")
#    os.system("echo 'TPROP = T' >> KPATH.phonopy ")
    
#########################################PHONON
#    atoms_ase = read("POSCAR")

# --------------------------
# 2. 构建 Phonopy 原胞 & 超胞
# --------------------------
    unitcell = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        scaled_positions=atoms.get_scaled_positions()
    )

    supercell_matrix = np.diag([2, 2, 2])
    phonon = Phonopy(unitcell, supercell_matrix)
    phonon.generate_displacements(distance=0.01)

    supercells = phonon.supercells_with_displacements
    print(f"生成 {len(supercells)} 个位移超胞")

# --------------------------
# 3. 计算力（ASE + MACE）
# --------------------------
    forces = []

    for i, sc in enumerate(supercells):
        print(f"计算超胞位移 {i+1}/{len(supercells)}")
        ase_scell = Atoms(
            symbols=sc.symbols,
            cell=sc.cell,
            scaled_positions=sc.scaled_positions,
            pbc=True
        )
        ase_scell.calc = calc
        f = ase_scell.get_forces()
        forces.append(f)

# --------------------------
# 4. 导入力并生成力常数
# --------------------------
# Phonopy 2.44.0 正确用法
# 直接生成力常数
    phonon.produce_force_constants(forces=forces)

# 保存 FORCE_SETS 文件
    dataset = phonon._dataset  # Phonopy 2.44.0 获取 dataset
    write_FORCE_SETS(dataset, filename="FORCE_SETS")
    print("生成力常数和 FORCE_SETS 成功")



# --------------------------
# 5. 声子能带计算
# --------------------------
# --------------------------
# 6. 声子 DOS
# --------------------------
    os.system("phonopy KPATH.phonopy")
    os.system("phonopy-bandplot --gnuplot band.yaml > band.dat")
# --------------------------
# 7. 热力学性质 F(T), S(T), Cv(T)
# --------------------------


# 运行振动分析
    vib = Vibrations(atoms, name='vib_mace', delta=delta)
    vib.run()

# 获取频率（单位 cm^-1）和模式能量（单位 eV）
    freqs = vib.get_frequencies()
    energies = vib.get_energies()  # = ħω（单位 eV）

# 去除零频（平移或旋转）
    mask = freqs > 1e-3
    freqs = freqs[mask]
    energies = energies[mask]

# 打印频率统计
    print(f"Total modes: {len(freqs)}")
    print(f"Lowest frequency: {freqs.min():.3f} cm^-1")
    print(f"Highest frequency: {freqs.max():.3f} cm^-1")

# 自由能计算函数
    kB = 8.617333262145e-5  # eV/K

    def F_vib(T):
        """Return vibrational free energy (eV) at temperature T (K)."""
        if T == 0:
            return 0.5 * np.sum(energies)  # 零点能
        beta = 1.0 / (kB * T)
        return np.sum(0.5 * energies + kB*T*np.log(1 - np.exp(-energies * beta)))

# 输出各温度自由能
    fileG = open(pos.split('.')[0]+'_G.txt','w')
    fileG.write("\nVibrational free energy (including ZPE):"+'\n')
    fileG.write(" T (K)    F_vib (eV)"+'\n')
    fileG.write("------------------------"+'\n')
    for T in temperatures:
        F = F_vib(T)
        fileG.write(f"{T:6.1f}   {F:10.6f}"+'\n')
        if T == T_ref:
            G = F.real + e
            
            fileRes = open(pos.split('.')[0]+'_G.res','w')
            fileRes.write('TITL '+pos+'  '+str(p)+'  '+str(v)+'  '+str(G)+' 0 0 '+str(TOT_Atoms)+' ('+SYS+')  n - 1'+'\n')
            fileRes.write('REM'+'\n')
            fileRes.write('REM'+'\n')
            fileRes.write('CELL 1.54180    '+str(round(a,6))+'    '+str(round(b,6))+'    '+str(round(c,6))+'    '+str(round(alpha,6))+'    '+str(round(beta,6))+'    '+str(round(gamma,6))+'    '+'\n')
            fileRes.write('LATT -1'+'\n')
            fileRes.write('SFAC ')
            for elem, num in counts.items():
                fileRes.write(elem+'  ')
            fileRes.write('\n')
            for j, atom in enumerate(atoms):
                x, y, z = atoms.get_scaled_positions()[j]
                fileRes.write(str(atom.symbol)+'  ')
                fileRes.write(str(1)+'  ')
                fileRes.write(str(x)+'  ')
                fileRes.write(str(y)+'  ')
                fileRes.write(str(z)+'  ')
                fileRes.write('1.0  ')
                fileRes.write('\n')
            fileRes.write('END'+'\n')
            fileRes.close()

    fileG.close()
    os.system('mv '+pos.split('.')[0]+'_G.res '+pos.split('.')[0]+'_G.txt' +' G')
# 可选：保存频率列表
#np.savetxt("vib_frequencies_cm-1.txt", freqs, header="Frequencies (cm^-1)")

# 清理临时文件（默认 Vibrations 会生成 vib.* 文件夹）
    os.system('rm -rf vib_mace*')





    os.system('mv mesh.yaml HIGH_SYMMETRY_POINTS KPATH.phonopy POSCAR FORCE_SETS force.yaml band.yaml band.dat total_dos.dat thermal_properties.yaml projected_dos.dat  phonon/'+pos.split('.')[0])





