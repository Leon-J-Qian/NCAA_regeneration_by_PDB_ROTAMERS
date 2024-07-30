import os
import argparse
import re
import numpy as np
import shutil
from rdkit.Chem import Lipinski
from Bio.PDB import PDBParser, PDBIO, Superimposer
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
import random

# 定义乙酰化和甲氨基化的SMARTS模式
acetylation_smarts = '[N:1][C:2][C:3](=[O:4])>>[C:5][C:7](=[O:6])[N:1][C:2][C:3](=[O:4])'
amidation_smarts = '[C:1][C:2](=[O:3])[O:4]>>[C:1][C:2]([N:5][C:6])(=[O:3])'

def count_rotatable_bonds(smiles_file, cut_off):
    with open(smiles_file, 'r') as file:
        smiles_list = [line.strip() for line in file.readlines()]

    results = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"SMILES字符串 '{smiles}' 无法解析")
            continue

        num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
        if num_rotatable_bonds==2 and cut_off==10000:
            v_value=10
            c_value=10000
        elif num_rotatable_bonds==3 and cut_off==10000:
            v_value=30
            c_value=10000
        elif num_rotatable_bonds==4 and cut_off==10000:
            v_value=90
            c_value=10000
        elif num_rotatable_bonds==5 and cut_off==10000:
            v_value=270
            c_value=10000
        elif num_rotatable_bonds>=6 and cut_off==10000:
            v_value=810
            c_value=10000
        else:
            v_value=cut_off
            c_value=10000
        results.append((smiles, num_rotatable_bonds, c_value, v_value))

    return results

def process_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Unable to parse SMILES: '{smiles}'")

    # 创建乙酰化和甲氨基化的反应
    rxn_acetylation = AllChem.ReactionFromSmarts(acetylation_smarts)
    rxn_amidation = AllChem.ReactionFromSmarts(amidation_smarts)

    # 应用甲氨基化反应
    products_amidation = rxn_amidation.RunReactants((mol,))
    if not products_amidation:
        raise ValueError("C-terminal amidation failed")

    product_amidation = products_amidation[0][0]

    # 应用乙酰化反应
    products_acetylation = rxn_acetylation.RunReactants((product_amidation,))
    if not products_acetylation:
        raise ValueError("N-terminal acetylation failed")

    product_acetylation = products_acetylation[0][0]

    # 获取乙酰化后的分子的 SMILES
    final_smiles = Chem.MolToSmiles(product_acetylation, canonical=False)
    return final_smiles

def smiles_to_pdb(smiles_list, names_list, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if len(names_list) < len(smiles_list):
        names_list.extend(['UAA'] * (len(smiles_list) - len(names_list)))
    elif len(names_list) > len(smiles_list):
        raise ValueError("The number of provided names exceeds the number of SMILES strings.")

    generated_files = []
    for i, (smiles, name) in enumerate(zip(smiles_list, names_list)):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"SMILES字符串 '{smiles}' 无法解析")
            continue
        mol = Chem.AddHs(mol)
        # 设置随机种子
        random_seed = random.randint(1, 1000000)
        params = AllChem.ETKDG()
        params.randomSeed = random_seed
        AllChem.EmbedMolecule(mol, params)
        
        # 使用UFF力场进行结构优化
        ff = UFFGetMoleculeForceField(mol)
        ff.Initialize()
        ff.Minimize()

        # 获取结构优化后的能量值
        energy = ff.CalcEnergy()

        # 生成PDB文件名
        pdb_filename = os.path.join(output_dir, f'{name}_{i+1}.pdb')
        
        # 写入PDB文件并添加注释
        with open(pdb_filename, 'w') as pdb_file:
            pdb_file.write(f"REMARK SMILES: {smiles}\n")
            pdb_file.write(f"REMARK Energy after optimization: {energy:.2f} kcal/mol\n")
            pdb_file.write(Chem.rdmolfiles.MolToPDBBlock(mol))
        
        generated_files.append(pdb_filename)
        print(f"生成 {pdb_filename}")
        

    return generated_files

def gen_smiles(input_file, names_list):
    # Read input SMILES file and process them
    with open(input_file, 'r') as f:
        smiles_list = [line.strip() for line in f.readlines()]

    if len(names_list) != len(smiles_list):
        raise ValueError("The number of provided names does not match the number of SMILES strings.")

    processed_smiles_list = []
    for smiles in smiles_list:
        try:

            processed_smiles = process_smiles(smiles)
            processed_smiles_list.append(processed_smiles)
        except ValueError as e:
            print(e)
            continue
    return processed_smiles_list

def gen_conformers(processed_smiles_list, names_list):

    # Convert processed SMILES to PDB files
    output_dir = 'pdb_files'
    generated_files = smiles_to_pdb(processed_smiles_list, names_list, output_dir)

def append_pdb_to_combined_file(pdb_file, combined_file):
    with open(pdb_file, 'r') as f_pdb:
        pdb_content = f_pdb.read()
    
    with open(combined_file, 'a') as f_combined:
        f_combined.write(pdb_content)

def run_repeatedly(input_file, names_list):
    combined_file = 'combined_pdb_files.pdb'
    
    if os.path.exists(combined_file):
        os.remove(combined_file)

    processed_smiles_list = gen_smiles(input_file, names_list)

    for i in range(c_value):
        print(f"Running iteration {i + 1}...")
        gen_conformers(processed_smiles_list, names_list)
        
        # 将pdb_files文件夹中的所有PDB文件内容追加到combined_pdb_files.pdb中
        for pdb_file in os.listdir('pdb_files'):
            if pdb_file.endswith('.pdb'):
                pdb_filepath = os.path.join('pdb_files', pdb_file)
                append_pdb_to_combined_file(pdb_filepath, combined_file)
        
        print(f"Iteration {i + 1} completed.")

def extract_energy_from_pdb(pdb_content):
    lines = pdb_content.splitlines()
    energy = None
    for line in lines:
        if line.startswith("REMARK Energy after optimization:"):
            energy = float(re.search(r'\d+\.\d+', line).group())
            break
    if energy is None:
        raise ValueError("No energy information found in PDB content")
    return energy

def process_pdb_file(filepath, n_value, output_filepath, energy_value):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # 添加能量值注释行
    energy_line = f"REMARK Energy after optimization: {energy_value:.2f}\n"

    # 过滤掉不是以 "HETATM" 或 "ATOM  " 开头的行
    lines = [line for line in lines if line.startswith("HETATM") or line.startswith("ATOM  ")]

    # 字典存储对应的行
    line_dict = {"C2": None, "O1": None, "C1": None, "H1": None, "H2": None, "H3": None,
                 "N2": None, "C5": None, "H6": None, "H7": None, "H8": None, "H9": None,
                 "N1": None, "C3": None, "C4": None, "O2": None, "H4": None, "H5": None}
    other_lines = []

    for line in lines:
        atom_name = line[12:16].strip()
        if atom_name in line_dict:
            line_dict[atom_name] = line
        else:
            other_lines.append(line)

    # 按指定顺序重新排列行
    reordered_lines = []
    for key in ["C2", "O1", "C1", "H1", "H2", "H3", "N2", "C5", "H6", "H7", "H8", "H9", "N1", "C3", "C4", "O2"]:
        if line_dict[key] is not None:
            reordered_lines.append(line_dict[key])

    reordered_lines.extend(other_lines)

    if line_dict["H4"] is not None:
        reordered_lines.append(line_dict["H4"])

    if line_dict["H5"] is not None:
        reordered_lines.append(line_dict["H5"])

    # 修改第四列
    for i in range(len(reordered_lines)):
        if i < 6:
            reordered_lines[i] = reordered_lines[i][:17] + "ACE" + reordered_lines[i][20:]
        elif i < 12:
            reordered_lines[i] = reordered_lines[i][:17] + "NME" + reordered_lines[i][20:]
        else:
            reordered_lines[i] = reordered_lines[i][:17] + n_value + reordered_lines[i][20:]

    # 在开头添加能量值注释行
    reordered_lines.insert(0, energy_line)

    # 添加 "END" 行
    reordered_lines.append("END\n")

    # 写回新文件
    with open(output_filepath, 'w') as file:
        file.writelines(reordered_lines)

def split_combined_pdb(combined_pdb_file, output_dir, n_value):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pdb_files = []
    energies = {}

    with open(combined_pdb_file, 'r') as f:
        content = f.read().strip().split('END\n')

    for idx, pdb_block in enumerate(content):
        if pdb_block.strip():
            pdb_content = pdb_block.strip() + '\nEND'
            energy = extract_energy_from_pdb(pdb_content)

            # 创建PDB文件名，使用能量值命名
            pdb_filename = os.path.join(output_dir, f'conformer_{energy:.2f}.pdb')

            with open(pdb_filename, 'w') as pdb_out:
                pdb_out.write(pdb_content)

            pdb_files.append((pdb_filename, energy))
            energies[pdb_filename] = energy

            # 对每个拆分的PDB文件进行内容上的修改
            output_filepath = f"{pdb_filename}_processed.pdb"
            process_pdb_file(pdb_filename, n_value, output_filepath, energy)
            print(f"已处理PDB文件: {output_filepath}")

    # 按能量排序生成的PDB文件
    sorted_files = sorted(pdb_files, key=lambda x: x[1])

    # 创建新的合并PDB文件并写入排序后的内容
    combined_sorted_file = 'combined_sorted_pdb_files.pdb'
    with open(combined_sorted_file, 'w') as f_combined:
        for pdb_file, _ in sorted_files:
            with open(pdb_file + "_processed.pdb", 'r') as f_pdb:
                pdb_content = f_pdb.read()
                f_combined.write(pdb_content)

    return combined_sorted_file

# 定义主链原子名
backbone_atoms = {'C3', 'C4', 'O1', 'N1', 'O2'}

def get_backbone_atoms(structure):
    """获取主链原子"""
    return [atom for atom in structure.get_atoms() if atom.get_name() in backbone_atoms]

def align_and_overwrite(ref_structure, pdb_file, parser):
    """对齐结构并覆盖原始PDB文件"""
    structure = parser.get_structure(pdb_file, pdb_file)
    
    # 获取主链原子
    ref_atoms = get_backbone_atoms(ref_structure)
    mobile_atoms = get_backbone_atoms(structure)
    
    # 使用 Superimposer 进行对齐
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, mobile_atoms)
    super_imposer.apply(structure.get_atoms())

    # 保存对齐后的结构到原始PDB文件中
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def kabsch_algorithm(P, Q):
    """Kabsch算法计算最佳旋转矩阵"""
    # 中心化
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    
    # 计算协方差矩阵
    H = np.dot(P_centered.T, Q_centered)
    
    # 计算SVD
    U, S, Vt = np.linalg.svd(H)
    
    # 计算旋转矩阵
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

def get_coordinates(structure, atoms_to_exclude):
    """从PDB文件中提取排除指定原子外的坐标和名称"""
    atoms = []
    coords = []
    for atom in structure.get_atoms():
        if atom.get_name() not in atoms_to_exclude and 'H' not in atom.get_name():  # 排除指定原子和所有氢原子
            atoms.append(atom.get_name())
            coords.append(atom.get_coord())
    
    return np.array(coords), atoms

def align_and_calculate_rmsd(ref_structure, temp_structure, original_pdb_file):
    """计算两个结构之间的RMSD"""
    # 指定要排除的原子名称和所有氢原子
    atoms_to_exclude = {'C1', 'C2', 'C4', 'C5', 'O1', 'O2', 'N1', 'N2'}
    
    coords1, atom_names1 = get_coordinates(ref_structure, atoms_to_exclude)
    coords2, atom_names2 = get_coordinates(temp_structure, atoms_to_exclude)
    
    if coords1.shape != coords2.shape:
        raise ValueError("两个PDB文件的坐标数量不一致")
    
    rmsd_value = calculate_rmsd(coords1, coords2)
    
    return rmsd_value

def align_structures_and_check_rmsd(pdb_files, parser, aligned_combined_pdb_filename, rmsd_threshold):
    aligned_structures = []
    conformer_count = 0

    with open(aligned_combined_pdb_filename, 'w') as aligned_combined_pdb:
        ref_pdb_file = pdb_files[0]
        ref_structure = parser.get_structure('reference', ref_pdb_file)
        align_and_overwrite(ref_structure, ref_pdb_file, parser)
        aligned_structures.append(ref_structure)

        for pdb_file in pdb_files[1:]:
            if conformer_count >=v_value or conformer_count >= len(pdb_files):
                break

            align_and_overwrite(ref_structure, pdb_file, parser)
            temp_structure = parser.get_structure('temp_aligned', pdb_file)

            # 比较temp_structure与aligned_combined_pdb_files.pdb中所有结构的RMSD值
            should_keep = True
            for aligned_structure, aligned_pdb_file in zip(aligned_structures, pdb_files):
                rmsd = align_and_calculate_rmsd(aligned_structure, temp_structure, aligned_pdb_file)
                print(f"RMSD between {os.path.basename(pdb_file)} and {os.path.basename(aligned_pdb_file)}: {rmsd}")
                if rmsd < rmsd_threshold:
                    should_keep = False
                    break
            
            # 如果不存在RMSD小于阈值的结果，则将temp_structure写入最终的合并文件
            if should_keep:
                with open(pdb_file, 'r') as temp_file:
                    lines = temp_file.readlines()
                    for line in lines:
                        if line.startswith(('REMARK', 'ATOM', 'HETATM', 'END')):
                            aligned_combined_pdb.write(line)
                
                conformer_count += 1
                aligned_structures.append(temp_structure)

    # 输出实际生成的构象个数
    print(f"实际生成的构象个数为: {conformer_count}")

def screen_rmsd(rmsd_threshold):
    parser = PDBParser(QUIET=True)

    # 读取combined_pdb_files.pdb文件
    with open('combined_sorted_pdb_files.pdb', 'r') as f:
        content = f.read().strip().split('END\n')

    # 将每个PDB结构写入单独的文件
    pdb_files = []
    for i, pdb_content in enumerate(content):
        if pdb_content.strip():
            pdb_filename = f'pdb_{i}.pdb'
            with open(pdb_filename, 'w') as pdb_file:
                pdb_file.write(pdb_content + 'END\n')
            pdb_files.append(pdb_filename)

    # 创建一个新的文件来存储所有对齐后的PDB结构
    aligned_combined_pdb_filename = 'aligned_combined_pdb_files.pdb'

    align_structures_and_check_rmsd(pdb_files, parser, aligned_combined_pdb_filename, rmsd_threshold)

    # 删除中间生成的PDB文件
    for pdb_file in pdb_files:
        os.remove(pdb_file)

    print(f"已创建合并PDB文件: aligned_combined_pdb_files.pdb")

def delete_files():
    files_to_delete = ['combined_pdb_files.pdb', 'combined_sorted_pdb_files.pdb', 'aligned_combined_pdb_files.pdb']
    
    for file in files_to_delete:
        try:
            os.remove(file)
            print(f"Successfully deleted {file}")
        except FileNotFoundError:
            print(f"File {file} not found")
        except PermissionError:
            print(f"Permission denied to delete {file}")
        except Exception as e:
            print(f"Error deleting {file}: {e}")

    if os.path.exists('pdb_files'):
        shutil.rmtree('pdb_files')
        print("Deleted 'pdb_files' directory")

    if os.path.exists('split_pdb_files'):
        shutil.rmtree('split_pdb_files')
        print("Deleted 'split_pdb_files' directory")

# Define the directory containing split PDB files
split_pdb_directory = 'split_pdbs'

def remove_ace_nme_lines(pdb_filename):
    """Remove lines containing 'ACE' or 'NME' from a PDB file and renumber the second column."""
    lines = []
    current_residue_number = 1

    with open(pdb_filename, 'r') as f:
        for line in f:
            if 'ACE' in line or 'NME' in line:
                continue
            if line.startswith('ATOM'):
                # Re-number the second column
                line = f"{line[:6]}{current_residue_number:>4}{line[10:]}"
                current_residue_number += 1
            lines.append(line)

    # Write back to the original PDB file
    with open(pdb_filename, 'w') as f:
        f.write(''.join(lines))

def split_pdb_file(input_pdb_file):
    """Split a PDB file into individual conformations."""
    output_directory = 'split_pdbs'
    os.makedirs(output_directory, exist_ok=True)

    conformation_index = 1
    pdb_content = []

    with open(input_pdb_file, 'r') as f:
        for line in f:
            if line.strip() == 'END':
                pdb_filename = os.path.join(output_directory, f'conformation_{conformation_index}.pdb')
                with open(pdb_filename, 'w') as pdb_file:
                    pdb_file.write('\n'.join(pdb_content) + '\nENDMDL\n')

                conformation_index += 1
                pdb_content = []
            else:
                pdb_content.append(line.strip())

    if pdb_content:
        pdb_filename = os.path.join(output_directory, f'conformation_{conformation_index}.pdb')
        with open(pdb_filename, 'w') as pdb_file:
            pdb_file.write('\n'.join(pdb_content) + '\nENDMDL\n')

def process_all_pdb_files(directory, non_natural_aa_codes):
    """Process all PDB files in the given directory."""
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            pdb_filepath = os.path.join(directory, filename)
            remove_ace_nme_lines(pdb_filepath)
            # Replace [12:16] of each line with non-natural amino acid codes
            lines = []
            with open(pdb_filepath, 'r') as f:
                for idx, line in enumerate(f):
                    if line.strip() == 'ENDMDL':
                        lines.append(line)
                        continue
                    # Replace [12:16] of each line with the non-natural amino acid codes
                    if idx < len(non_natural_aa_codes):
                        code = non_natural_aa_codes[idx]
                        line = line[:12] + code + line[16:]
                    lines.append(line)
            
            # Write back to the original PDB file with modified lines
            with open(pdb_filepath, 'w') as f:
                f.write(''.join(lines))

            print(f"Processed: {filename}")

def modify(non_natural_aa_abbreviation):
    # Determine the params filename based on the three-letter abbreviation
    params_filename = f"{non_natural_aa_abbreviation}_temps.params"

    # Read XXX_temps.params file
    non_natural_aa_codes = []
    with open(params_filename, 'r') as params_file:
        for line in params_file:
            if line.startswith('ATOM'):
                code = line.strip()[5:9]  # Extract the [5:9] substring
                non_natural_aa_codes.append(code)

    # Split the input PDB file into individual conformations
    input_pdb_file = 'aligned_combined_pdb_files.pdb'
    split_pdb_file(input_pdb_file)

    # Process all PDB files in the split_pdb_directory
    process_all_pdb_files(split_pdb_directory, non_natural_aa_codes)

    # Combine split PDB files into a single merged PDB file
    combined_pdb_filename = 'merged_combined_pdb_files.pdb'
    combine_pdb_files(combined_pdb_filename)

    # Process combined PDB file to adjust HETATM records
    adjust_hetatm_records(combined_pdb_filename)

    # Remove split_pdb_directory
    remove_split_pdb_directory(split_pdb_directory)
    print(f"Combined PDB files into '{combined_pdb_filename}' and deleted '{split_pdb_directory}'.")

def combine_pdb_files(output_filename):
    """Combine all PDB files in split_pdb_directory into a single output file."""
    with open(output_filename, 'w') as output_file:
        conformation_index = 1
        for filename in sorted(os.listdir(split_pdb_directory)):
            if filename.endswith('.pdb'):
                pdb_filepath = os.path.join(split_pdb_directory, filename)
                with open(pdb_filepath, 'r') as pdb_file:
                    output_file.write(f"MODEL     {conformation_index}\n")
                    for line in pdb_file:
                        output_file.write(line)
                    conformation_index += 1

def remove_split_pdb_directory(directory):
    """Remove the split_pdb_directory."""
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

    # Remove the directory itself
    try:
        os.rmdir(directory)
    except Exception as e:
        print(f"Failed to delete {directory}. Reason: {e}")

def adjust_hetatm_records(pdb_filename):
    """Adjust HETATM records by subtracting 12 from the second column (index 10:11) and adding the difference back to the original position."""
    lines = []
    with open(pdb_filename, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                original_number = int(line[9:11])
                modified_number = original_number - 12
                # Adjust right alignment
                if modified_number < 10:
                    line = f"{line[:9]} {modified_number}{line[11:]}"
                else:
                    line = f"{line[:9]}{modified_number}{line[11:]}"
            lines.append(line)

    with open(pdb_filename, 'w') as f:
        f.writelines(lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process SMILES and convert to PDB files.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input SMILES file.')
    parser.add_argument('-n', '--names', required=True, nargs='+', help='Names for each SMILES (UAA if not specified).')
    parser.add_argument('-m', '--params', required=True, help='Path to the params file.')
    parser.add_argument('-d', '--rmsd_threshold', type=float, default=0.001,
                        help='RMSD threshold for filtering conformers')
    parser.add_argument('-c', '--cut_off', type=int, default=10000,
                        help='Cut off value for the number of conformers')
    args = parser.parse_args()

    results = count_rotatable_bonds(args.input, args.cut_off)
    for smiles, num_rotatable_bonds, c_value, v_value in results:
        print(f"SMILES: {smiles}  自由旋转键数目: {num_rotatable_bonds}  c_value: {c_value}  v_value: {v_value}")
        # 将c_value和v_value输出为全局变量
        globals()['c_value'] = c_value
        globals()['v_value'] = v_value

    try:
        run_repeatedly(args.input, args.names)
    except Exception as e:
        print(f"Error: {e}")
        
    n_value = args.names[0]
    combined_pdb_file = 'combined_pdb_files.pdb'  # 假设这是您的合并PDB文件名
    output_directory = 'split_pdb_files'  # 存放拆分后PDB文件的目录
    try:
        sorted_combined_file = split_combined_pdb(combined_pdb_file, output_directory, n_value)
        print(f"已创建排序后的合并PDB文件: {sorted_combined_file}")

    except Exception as e:
        print(f"错误: {e}")
            
    rmsd_threshold = args.rmsd_threshold
    if rmsd_threshold==0.001 and v_value==10:
        rmsd_threshold = 0.80
    elif rmsd_threshold==0.001 and v_value==30:
        rmsd_threshold = 0.60
    elif rmsd_threshold==0.001 and v_value==90:
        rmsd_threshold = 0.40
    elif rmsd_threshold==0.001 and v_value==270:
        rmsd_threshold = 0.20
    elif rmsd_threshold==0.001 and v_value==810:
        rmsd_threshold = 0.10
    print(rmsd_threshold)
    try:
        screen_rmsd(rmsd_threshold)

    except Exception as e:
        print(f"错误: {e}")
    
    non_natural_aa_abbreviation = n_value

    try:
        modify(non_natural_aa_abbreviation)

    except Exception as e:
        print(f"Error: {e}")
    delete_files()
        