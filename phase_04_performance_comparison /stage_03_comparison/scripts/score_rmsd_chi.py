import os
import argparse
import numpy as np
from Bio import PDB
import pandas as pd

def find_lowest_energy_pdb(score_file):
    min_energy = float('inf')
    min_pdb = None

    with open(score_file, 'r') as file:
        lines = file.readlines()[2:]  # 跳过前两行
        for line in lines:
            columns = line.split()
            energy = float(columns[1])
            pdb_name = columns[-1]
            if energy < min_energy:
                min_energy = energy
                min_pdb = pdb_name
    
    return min_pdb, min_energy

def kabsch_algorithm(P, Q):
    """Kabsch算法计算最佳旋转矩阵"""
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    H = np.dot(P_centered.T, Q_centered)
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    d = np.sign(np.linalg.det(np.dot(V, U.T)))
    R = np.dot(V, np.dot(np.diag([1, 1, d]), U.T))
    return R

def calculate_rmsd(P, Q):
    """计算两个点集之间的RMSD"""
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    R = kabsch_algorithm(P_centered, Q_centered)
    P_rot = np.dot(P_centered, R)
    diff = P_rot - Q_centered
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

def get_heavy_atom_coordinates(pdb_file, residue_name):
    """从PDB文件中提取指定残基的所有重原子的坐标"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    heavy_atom_coordinates = []
    for residue in structure.get_residues():
        if residue.get_resname() == residue_name:
            for atom in residue.get_atoms():
                if atom.element != 'H':  # 排除氢原子
                    heavy_atom_coordinates.append(atom.get_coord())
    return np.array(heavy_atom_coordinates)

def split_conformers(pdb_file):
    """将包含多个构象的PDB文件拆分成单个构象的PDB文件"""
    with open(pdb_file, 'r') as file:
        content = file.read()
    conformers = content.split('END\n')
    pdb_files = []
    for i, conformer in enumerate(conformers):
        if conformer.strip():
            temp_pdb = f"temp_conformer_{i+1}.pdb"
            with open(temp_pdb, 'w') as file:
                file.write(conformer + 'END\n')
            pdb_files.append(temp_pdb)
    return pdb_files

def calculate_rmsd_for_pdbs(single_pdb, conformers_pdb, residue_name):
    """计算两个PDB文件的RMSD，其中一个包含一个构象，另一个包含多个构象"""
    single_coords = get_heavy_atom_coordinates(single_pdb, residue_name)
    temp_pdb_files = split_conformers(conformers_pdb)
    rmsd_list = []
    for i, temp_pdb in enumerate(temp_pdb_files):
        conformer_coords = get_heavy_atom_coordinates(temp_pdb, residue_name)
        if single_coords.shape[0] != conformer_coords.shape[0]:
            raise ValueError(f"原子数量不一致: {single_pdb} 有 {single_coords.shape[0]} 个重原子，{temp_pdb} 有 {conformer_coords.shape[0]} 个重原子")
        rmsd = calculate_rmsd(single_coords, conformer_coords)
        rmsd_list.append((i+1, rmsd))
        os.remove(temp_pdb)
    return rmsd_list

def merge_pdb_files(n, i):
    """合并PDB文件"""
    output_file = "merged.pdb"
    with open(output_file, 'w') as merged_file:
        for filename in os.listdir():
            if filename.startswith(f"{i}_opt") and filename.endswith(".pdb"):
                with open(filename, 'r') as pdb_file:
                    lines = pdb_file.readlines()
                    for line in lines:
                        if line[17:20].strip() == n:
                            merged_file.write(line)
                    merged_file.write("END\n")
    return output_file

def read_chi_file(chi_file):
    chi_values = {}
    with open(chi_file, 'r') as file:
        for line in file:
            key, value = line.strip().split('=')
            chi_values[key.strip()] = float(value.strip())
    return chi_values

def calculate_dihedral(coords):
    p0, p1, p2, p3 = coords
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    dihedral_angle = np.degrees(np.arctan2(y, x))
    return dihedral_angle
    
def read_params_file(params_file):
    chi_list = []
    with open(params_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('CHI') and 'H' not in line[3:]:
                chi_list.append(line.strip())
    return chi_list

def process_chi_list(chi_list):
    res_chi = {}
    for i, chi in enumerate(chi_list):
        elements = chi.split()
        res_chi['res_chi{}'.format(i+1)] = ' '.join(elements[2:6])
    return res_chi

def process_pdb_files(res_chi, amino_acid, chi_values):
    results = []
    for file in os.listdir('.'):
        if '_opt_' in file and file.endswith('.pdb'):
            chi_values_calculated = []
            with open(file, 'r') as pdb:
                lines = pdb.readlines()
                for key, res in res_chi.items():
                    coords = []
                    atoms = res.split()
                    for atom in atoms:
                        for line in lines:
                            if (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20].strip() == amino_acid and line[12:16].strip() == atom:
                                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                                coords.append(np.array([x, y, z]))
                                break
                    if len(coords) == 4:
                        chi_angle = calculate_dihedral(coords)
                        chi_values_calculated.append(chi_angle)
            deviations = []
            for i, chi_angle in enumerate(chi_values_calculated):
                deviation = abs(chi_angle - chi_values['chi{}'.format(i+1)])  # Absolute deviation
                deviations.append('{:.3f}'.format(deviation))  # Format to 3 decimal places
            results.append([file] + deviations)
    return results

def write_results(results):
    with open('CHI_results.csv', 'w') as file:
        for result in results:
            file.write(', '.join(result) + '\n')

def main():
    parser = argparse.ArgumentParser(description="合并PDB文件并计算RMSD及CHI角偏差")
    parser.add_argument("-n", required=True, help="要筛选的残基名称")
    parser.add_argument("-i", required=True, help="要筛选的文件前缀")
    parser.add_argument("-r", required=True, help="参考的PDB文件")
    parser.add_argument("-s", default="score_opt.sc", help="score文件路径")
    parser.add_argument('-c', default="chi.txt", help='Input chi.txt file')
    args = parser.parse_args()

    # 找到能量最低的PDB文件
    pdb_name, energy = find_lowest_energy_pdb(args.s)
    print(f'The PDB with the lowest energy is {pdb_name} with an energy of {energy}')

    # 合并PDB文件
    merged_pdb = merge_pdb_files(args.n, args.i)

    # 计算所有PDB文件的RMSD
    rmsd_results = calculate_rmsd_for_pdbs(args.r, merged_pdb, args.n)
    
    # 保存到CSV文件
    df = pd.DataFrame(rmsd_results, columns=['Conformer', 'RMSD'])
    df.to_csv('rmsd_results.csv', index=False)
    
    # 计算平均值和标准差
    mean_rmsd = df['RMSD'].mean()
    std_rmsd = df['RMSD'].std()
    
    # 输出平均值和标准差
    print(f"RMSD平均值: {mean_rmsd:.3f}")
    print(f"RMSD标准差: {std_rmsd:.3f}")

    # 计算能量最低的PDB文件的RMSD
    lowest_pdb_file = f"{pdb_name}.pdb"
    lowest_rmsd_results = calculate_rmsd_for_pdbs(args.r, lowest_pdb_file, args.n)
    
    # 找到能量最低的PDB的RMSD值
    lowest_rmsd = lowest_rmsd_results[0][1]  # 因为只有一个构象
    print(f'The RMSD of the lowest energy PDB ({pdb_name}) is {lowest_rmsd:.3f}')
    
    # 处理chi角
    chi_values = read_chi_file(args.c)
    params_file = f'{args.n}.params'
    params_lines = read_params_file(params_file)
    res_chi = process_chi_list(params_lines)
    chi_results = process_pdb_files(res_chi, args.n, chi_values)
    write_results(chi_results)
    
    # 输出chi角偏差
    for result in chi_results:
        print(f'PDB File: {result[0]}, ' + ', '.join([f'Chi{i+1} Deviation: {result[i+1]}' for i in range(len(result) - 1)]))
        
    # 输出能量最低的PDB文件的chi角偏差
    for result in chi_results:
        if result[0] == f"{pdb_name}.pdb":
            print(f'能量最低的PDB文件 ({pdb_name}) 的CHI角偏差:')
            print(', '.join([f'Chi{i+1}偏差: {result[i+1]}' for i in range(len(result) - 1)]))
            break
    
    # 写入结果到文件
    with open('results.txt', 'w') as file:
        file.write(f'The PDB with the lowest energy is {pdb_name} with an energy of {energy}\n')
        file.write(f"RMSD平均值: {mean_rmsd:.3f}\n")
        file.write(f"RMSD标准差: {std_rmsd:.3f}\n")
        file.write(f'The RMSD of the lowest energy PDB ({pdb_name}) is {lowest_rmsd:.3f}\n')
        file.write(f'能量最低的PDB文件 ({pdb_name}) 的CHI角偏差:\n')
        for result in chi_results:
            if result[0] == f"{pdb_name}.pdb":
                file.write(', '.join([f'Chi{i+1}偏差: {result[i+1]}' for i in range(len(result) - 1)]) + '\n')
                break
                
if __name__ == "__main__":
    main()