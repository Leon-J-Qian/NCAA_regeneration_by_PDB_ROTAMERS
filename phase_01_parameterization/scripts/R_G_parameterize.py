import os
import argparse
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg
import time
import shutil

# 定义乙酰化和甲氨基化的SMARTS模式
acetylation_smarts = '[N:1][C:2][C:3](=[O:4])>>[C:5][C:7](=[O:6])[N:1][C:2][C:3](=[O:4])'    # 乙酰化，为主链N连上乙酰基[C:5][C:7](=[O:6])
amidation_smarts = '[C:1][C:2](=[O:3])[O:4]>>[C:1][C:2]([N:5][C:6])(=[O:3])'                 # 甲氨基化，为羧基主链C连上甲氨基([N:5][C:6])(=[O:3]。注意这里reactant里需要将单键O[O:4]标注出来，不然会报错

#根据smiles的电荷信息计算体系净电荷值
def calculate_system_charge(smiles_filepath):
    with open(smiles_filepath, 'r') as file:
        smiles = file.read().strip()
    
    negative_count = smiles.count('-')
    positive_count = smiles.count('+')
    
    system_charge = positive_count - negative_count
    return system_charge
    
def process_smiles(smiles):
    # 识别输入的smiles是否可解析
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

    # 获取封端后的分子的 SMILES
    final_smiles = Chem.MolToSmiles(product_acetylation, canonical=False)
    print(final_smiles)
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
        AllChem.EmbedMolecule(mol)
        
        # 生成PDB文件名
        pdb_filename = os.path.join(output_dir, f'{name}_{i+1}.pdb')
        
        # 写入PDB文件并添加注释(该NCAA封端后的smiles，便于使用者检查)
        with open(pdb_filename, 'w') as pdb_file:
            pdb_file.write(f"REMARK {smiles}\n")
            pdb_file.write(Chem.rdmolfiles.MolToPDBBlock(mol))
        
        generated_files.append(pdb_filename)
        print(f"生成 {pdb_filename}")

        # 读取并打印PDB文件的内容
        with open(pdb_filename, 'r') as pdb_file_print:
            pdb_content = pdb_file_print.read()
            print(pdb_content)

    return generated_files

def process_pdb_file(filepath, n_value, output_filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # 过滤掉不是以 "HETATM" 或 "ATOM  " 开头的行
    lines = [line for line in lines if line.startswith("HETATM") or line.startswith("ATOM  ")]

    # 字典存储需要进行调整的行
    line_dict = { "C2": None, "O1": None, "C1": None, "H1": None, "H2": None, "H3": None,
                  "N2": None, "C5": None, "H6": None, "H7": None, "H8": None, "H9": None,
                  "N1": None, "C3": None, "C4": None, "O2": None, "H4": None, "H5": None }
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

    #后置主链H原子
    if line_dict["H4"] is not None:
        reordered_lines.append(line_dict["H4"])

    if line_dict["H5"] is not None:
        reordered_lines.append(line_dict["H5"])

    # 修改氨基酸名称，指认出ACE和NME的部分
    for i in range(len(reordered_lines)):
        if i < 6:
            reordered_lines[i] = reordered_lines[i][:17] + "ACE" + reordered_lines[i][20:]
        elif i < 12:
            reordered_lines[i] = reordered_lines[i][:17] + "NME" + reordered_lines[i][20:]
        else:
            reordered_lines[i] = reordered_lines[i][:17] + n_value + reordered_lines[i][20:]

    # 添加 "END" 行
    reordered_lines.append("END\n")

    # 写回新文件
    with open(output_filepath, 'w') as file:
        file.writelines(reordered_lines)

# 遍历所有pdb，对其执行重排操作
def process_files(filepaths, n_values, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for i, (filepath, n_value) in enumerate(zip(filepaths, n_values)):
        output_filepath = os.path.join(output_dir, f'{n_value}_{i+1}.pdb')
        process_pdb_file(filepath, n_value, output_filepath)
        print(f"Processed {filepath}")

        # 读取并打印PDB文件的内容
        with open(output_filepath, 'r') as output_print:
            output_content = output_print.read()
            print(output_content)

#定义二面角设置函数，应用于adjust_dihedrals_in_pdb_files函数的二面角调整，在发生错误时提供错误信息(如果该步出现报错则很有可能是上一步原子重排没有排好)
def set_dihedral_angle(mol, atom_indices, angle):
    try:
        SetDihedralDeg(mol.GetConformer(), *atom_indices, angle)
    except Exception as e:
        print(f"Error setting dihedral angle for atoms {atom_indices} to {angle} degrees: {e}")
        
def adjust_dihedrals_in_pdb_files(filepaths):
    for pdb_path in filepaths:
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol is None:
            print(f"Unable to load PDB file: {pdb_path}")
            continue
        
        # 获取原子索引
        atom_indices = {atom.GetPDBResidueInfo().GetName().strip(): atom.GetIdx() for atom in mol.GetAtoms()}

        # 检查 N1 原子连接的非氢原子数目
        N1_atom = mol.GetAtomWithIdx(atom_indices['N1'])
        non_h_neighbors = [nbr for nbr in N1_atom.GetNeighbors() if nbr.GetSymbol() != 'H']
        
        #根据N原子连接的非H原子数判断对象是peptide还是peptoid，并根据情况设置其对应的优势二面角
        if len(non_h_neighbors) == 2:
            # peptide
            set_dihedral_angle(mol, [atom_indices['C2'], atom_indices['N1'], atom_indices['C3'], atom_indices['C4']], -150.0)
            set_dihedral_angle(mol, [atom_indices['N1'], atom_indices['C3'], atom_indices['C4'], atom_indices['N2']], 150.0)
        elif len(non_h_neighbors) == 3:
            # peptoid
            set_dihedral_angle(mol, [atom_indices['C2'], atom_indices['N1'], atom_indices['C3'], atom_indices['C4']], -120.0)
            set_dihedral_angle(mol, [atom_indices['N1'], atom_indices['C3'], atom_indices['C4'], atom_indices['N2']], 90.0)
        else:
            print(f"Unexpected number of non-hydrogen neighbors for N1 in {pdb_path}")
            continue
        
        # 写回调整后的PDB文件
        with open(pdb_path, 'w') as f:
            f.write(Chem.MolToPDBBlock(mol))

        print(f"Adjusted dihedrals in {pdb_path}")

#过滤掉CONECT行，键连信息的存在有时会导致antechamber转化gjf文件时报错
def remove_conect_lines_from_pdb(pdb_filepath):
    with open(pdb_filepath, 'r') as f:
        lines = f.readlines()

    # 过滤掉以 "CONECT" 开头的行
    lines = [line for line in lines if not line.startswith("CONECT")]

    # 写回新文件
    with open(pdb_filepath, 'w') as f:
        f.writelines(lines)

    print(f"Removed CONECT lines from {pdb_filepath}")

#提交gaussian作业
def submit_gaussian_job(gjf_filepath):
    command = f"/home/rotations/Leon/g16/g16/g16 < {gjf_filepath} > {gjf_filepath.replace('.gjf', '.log')} &"
    subprocess.Popen(command, shell=True)
    print(f"Submitted Gaussian job for {gjf_filepath}")
    
def modify_gjf_file(gjf_filepath, system_charge):
    with open(gjf_filepath, 'r') as file:
        lines = file.readlines()

    # 检索并修改电荷行，为gjf文件提供正确的电荷信息
    for i in range(len(lines)):
        if lines[i].strip() == "0   1":
            lines[i] = f"{system_charge}   1\n"
            break

    # 删除第一、第二行，并定义内存和核数。可根据需要自行更改
    if len(lines) > 1:
        lines[0] = "%mem=100GB\n"
        lines[1] = "%nprocshared=24\n"
    else:
        lines.insert(0, "%mem=100GB\n")
        lines.insert(1, "%nprocshared=24\n")

    # 在最后一行非空行后添加固定二面角信息，以确保结构优化时Phi，Psi角固定不动
    non_empty_lines = [line for line in lines if line.strip()]
    last_non_empty_line_index = lines.index(non_empty_lines[-1])
    lines.insert(last_non_empty_line_index + 1, '\n13 14 15 7 F\n1 13 14 15 F\n\n\n\n')

    with open(gjf_filepath, 'w') as file:
        file.writelines(lines)

    print(f"Modified {gjf_filepath}")

#生成Gaussian输入文件
def generate_gaussian_input(filepaths, output_dir, system_charge):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for pdb_path in filepaths:
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        gjf_filename = os.path.join(output_dir, f'{base_name}.gjf')
        command = f'antechamber -i {pdb_path} -fi pdb -o {gjf_filename} -fo gcrt -gk "# opt=(modredundant,loose) b3lyp/6-31+g(d) scrf=(solvent=water) empiricaldispersion=gd3bj"'
        os.system(command)

        modify_gjf_file(gjf_filename, system_charge)

        # 使用 sed 命令替换 Cl 为 C，以避免因为原子编号问题导致的报错。如果体系存在Cl原子则需要注释掉以下命令并手动检查原子编号是否正常
        sed_command = f'sed -i "s/Cl/ C/g" {gjf_filename}'
        os.system(sed_command)

        print(f"Generated {gjf_filename}")

        # 复制生成的.gjf文件为同名的.log文件
        log_filename = os.path.join(output_dir, f'{base_name}.log')
        shutil.copy(gjf_filename, log_filename)

        # 提交Gaussian作业
        submit_gaussian_job(gjf_filename)

#检查函数，用以监督Gaussian作业是否结束运行
def is_gaussian_job_completed(log_filepath):
    with open(log_filepath, 'r') as file:
        lines = file.readlines()
    return any("Normal termination" in line for line in lines)
    
def process_log_files(input_file, output_dir):
    with open(input_file, 'r') as f:
        smiles_list = [line.strip() for line in f.readlines()]
    
    log_files = [f for f in os.listdir(output_dir) if f.endswith('.log')]
    chirality_dict = {}

    for log_file in log_files:
        base_name = os.path.splitext(log_file)[0]
        idx = int(base_name.split('_')[-1]) - 1  # assuming the log files are named like 'name_1.log', 'name_2.log', etc.
        smiles = smiles_list[idx]

        #判断结构优化是否正常结束，以防止因非正常原因停止的Gaussian结构优化所输出的log文件被应用于后续的处理过程中
        if not is_gaussian_job_completed(os.path.join(output_dir, log_file)):
            print(f"Gaussian job for {log_file} is not yet completed. Skipping for now.")
            continue

        #判断氨基酸为L型还是D型
        if '[C@@H]' in smiles:
            chirality_dict[log_file] = 'L'
        elif '[C@H]' in smiles:
            chirality_dict[log_file] = 'D'
        else:
            print(f"Cannot determine chirality for {log_file} from SMILES: {smiles}")
            continue
        
        mol_output_dir = 'mol'
        if not os.path.exists(mol_output_dir):
            os.makedirs(mol_output_dir)
        
        #将log文件转化为mol文件
        log_filepath = os.path.join(output_dir, log_file)
        mol_filepath = os.path.join(mol_output_dir, f'{base_name}_opt.mol')
        command = f'obabel -i g16 {log_filepath} -o mol -O {mol_filepath}'
        os.system(command)
        
        with open(mol_filepath, 'r') as mol_file:
            mol_lines = mol_file.readlines()
        
        # 删除最后一行"M  END"
        if mol_lines[-1].strip() == "M  END":
            mol_lines = mol_lines[:-1]

        # 根据氨基酸类型添加原子指认信息
        if chirality_dict[log_file] == 'L':
            mol_lines.extend([
                "M  ROOT 13\n",
                "M  POLY_N_BB 13\n",
                "M  POLY_CA_BB 14\n",
                "M  POLY_C_BB 15\n",
                "M  POLY_O_BB 16\n",
                "M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12\n",
                "M  POLY_UPPER 7\n",
                "M  POLY_LOWER 1\n",
                "M  POLY_PROPERTIES PROTEIN L_AA ALPHA_AA\n",
                "M  END\n"
            ])
        elif chirality_dict[log_file] == 'D':
            mol_lines.extend([
                "M  ROOT 13\n",
                "M  POLY_N_BB 13\n",
                "M  POLY_CA_BB 14\n",
                "M  POLY_C_BB 15\n",
                "M  POLY_O_BB 16\n",
                "M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12\n",
                "M  POLY_UPPER 7\n",
                "M  POLY_LOWER 1\n",
                "M  POLY_PROPERTIES PROTEIN D_AA ALPHA_AA\n",
                "M  END\n"
            ])
        
        with open(mol_filepath, 'w') as mol_file:
            mol_file.writelines(mol_lines)

        print(f"Processed and saved {mol_filepath}")

        # 读取并打印mol文件的内容
        with open(mol_filepath, 'r') as mol_print:
            mol_content = mol_print.read()
            print(mol_content)

def molfile_to_params(mol_filepath, name):
    # 创建params文件的完整路径
    base_name = os.path.splitext(os.path.basename(mol_filepath))[0]
    params_filename = f"{base_name}.params"
    params_filepath = os.path.join(os.getcwd(), params_filename)

    # 检查参数文件是否已经存在，如果存在则skip
    if os.path.exists(params_filepath):
        print(f"Params file {params_filepath} already exists. Skipping conversion.")
        return
    
    # 通过命令行调用molfile_to_params_polymer.py脚本，执行参数化
    # 实际使用中需要根据当前操作环境下的脚本路径对以下命令进行修改
    command = f"python2 /home/fzwang/rosetta_bin_linux_2020.25.61318_bundle/main/demos/public/using_ncaas_protein_peptide_interface_design/HowToMakeResidueTypeParamFiles/scripts/molfile_to_params_polymer.py -n {name} --polymer {mol_filepath}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # 检查任务是否正确运行
    if result.returncode != 0:
        print(f"Error converting {mol_filepath} to params file. Command output:\n{result.stdout}\n{result.stderr}")
        return

    # 打印成功信息
    print(f"Converted {mol_filepath} to {params_filepath}")
    
def molfile_to_params_temps(mol_filepath, name):
    # 创建temps.params路径
    base_name = os.path.splitext(os.path.basename(mol_filepath))[0]
    params_filename = f"{base_name}_temps.params"
    params_filepath = os.path.join(os.getcwd(), params_filename)

    # 检查参数文件是否已经存在
    if os.path.exists(params_filepath):
        print(f"Params file {params_filepath} already exists. Skipping conversion.")
        return
    
    # 使用molfile_to_params_polymer_modify.py脚本进行参数化
    command = f"python2 /home/fzwang/rosetta_bin_linux_2020.25.61318_bundle/main/demos/public/using_ncaas_protein_peptide_interface_design/HowToMakeResidueTypeParamFiles/scripts/molfile_to_params_polymer_modify.py -n {name}_temps --no_reorder --polymer {mol_filepath}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # 检查任务是否被顺利运行
    if result.returncode != 0:
        print(f"Error converting {mol_filepath} to params file. Command output:\n{result.stdout}\n{result.stderr}")
        return

    # 打印运行成功信息
    print(f"Converted {mol_filepath} to {params_filepath}")
    
# 创建resp拟合的gjf文件并提交，运行完毕后使用ff14SB力场将其转化为mol2文件
def generate_opt(res, resp):

    os.system(f'antechamber -i {res}.log -fi gout -o {resp}.gjf -fo gcrt -gk "# HF/6-31G*  SCF=Tight  Pop=MK  iop(6/33=2,  6/41=10, 6/42=15)"')
    print(f'\nGaussian RESP-calculation input file ({resp}.gjf) for {res}.mol2 has already been generated by antechamber!!\n')

    os.system(f'/home/rotations/Leon/g16/g16/g16 {resp}.gjf && antechamber -i {resp}.log -fi gout -o {resp}.mol2 -fo mol2 -at amber -pf y -c resp')
    print(f'\nThe optimized structure with RESP charge has been output to {resp}_atomtype.mol2 and needs to be further processed!!\n')

def atom_type_adjust(resp_folder):
    # 遍历RESP文件夹中的mol2文件
    for mol2_file in os.listdir(resp_folder):
        if mol2_file.endswith('.mol2'):
            # 提取文件名前三个字符作为res对象
            res = mol2_file[:3]

            # 查找对应的params文件
            params_file = f'{res}_temps.params'
            if not os.path.exists(params_file):
                print(f"Warning: {params_file} not found for {mol2_file}")
                continue
            
            # 读取params文件中ATOM开头行，将其全部保存进一个列表
            params_atom_lines = []
            with open(params_file, 'r') as f_params:
                for line in f_params:
                    if line.startswith('ATOM'):
                        params_atom_lines.append(line)
            
            # 读取mol2文件内容
            with open(os.path.join(resp_folder, mol2_file), 'r') as f_mol2:
                mol2_lines = f_mol2.readlines()
            
            start_idx = -1
            end_idx = -1
            
            # 找到@<TRIPOS>ATOM和@<TRIPOS>BOND之间的行的索引范围
            try:
                start_idx = mol2_lines.index('@<TRIPOS>ATOM\n') + 1
                end_idx = mol2_lines.index('@<TRIPOS>BOND\n')
            except ValueError:
                print(f"Error: Unable to find '@<TRIPOS>ATOM' or '@<TRIPOS>BOND' in {mol2_file}")
                continue
    
            # 检索这些行中第一列>=13的行，并将它们保存进另一个列表
            mol2_atom_lines_to_replace = []
            for i in range(start_idx, end_idx):
                parts = mol2_lines[i].split()
                if len(parts) >= 1:
                    try:
                        atom_index = int(parts[0])
                    except ValueError:
                        continue  # 如果无法转换为整数，跳过该行
            
                    if atom_index >= 13:
                        mol2_atom_lines_to_replace.append((i, mol2_lines[i]))
    
            # 将params_atom_lines中的[5:9]部分覆盖mol2_atom_lines_to_replace的[7:11]部分。即使用params的原子名称来替代mol2文件中的原子名称
            for j, (i, line) in enumerate(mol2_atom_lines_to_replace):
                parts = list(line)
                if len(parts) >= 11 and j < len(params_atom_lines):
                    new_value = params_atom_lines[j][5:9]
                    parts[7:11] = new_value
                    mol2_lines[i] = ''.join(parts)
            
            # 写入更新后的mol2文件
            output_file_path = os.path.join(resp_folder, mol2_file)
            with open(output_file_path, 'w') as f_mol2:
                f_mol2.writelines(mol2_lines)
                
                    # 读取更新后的mol2文件内容
            with open(output_file_path, 'r') as f_updated_mol2:
                updated_mol2_lines = f_updated_mol2.readlines()
    #        print('\n\n\n\n\n\n\n',updated_mol2_lines)
            
            # 查找@<TRIPOS>ATOM和@<TRIPOS>BOND之间的行的索引范围
            try:
                start_idx = updated_mol2_lines.index('@<TRIPOS>ATOM\n') + 1
                end_idx = updated_mol2_lines.index('@<TRIPOS>BOND\n')
            except ValueError:
                print(f"Error: Unable to find '@<TRIPOS>ATOM' or '@<TRIPOS>BOND' in updated {mol2_file}")
                continue
            
            # 检查并调整列表中[7:8]为数字的元素，将数字移至字母后面，使原子名称符合gromacs字母在前数字在后的规范
            adjusted_lines = []
            for i in range(start_idx, end_idx):
                line = updated_mol2_lines[i]
                parts = list(line)
                if len(parts) >= 8 and parts[7].isdigit():
                    if len(parts) >= 11 and parts[9] != ' ' and parts[10] != ' ':
                        parts[11] = parts[7]
    #                    print(parts[11])
                    if len(parts) >= 10 and parts[9] != ' ':
                        parts[10] = parts[7]
    #                    print(parts[10])
                    if len(parts) >= 9 and parts[9] == ' ':
                        parts[9] = parts[7]
    #                    print(parts[9])
                    parts[7:8] = ' '
                    adjusted_lines.append(''.join(parts))
                else:
                    adjusted_lines.append(line)
    #                print('no adjust')
    #        print(updated_mol2_lines)
    #        print('\n\n\n\n\n\n\n',adjusted_lines)
            for i in range(start_idx, end_idx):
                mol2_lines[i] = adjusted_lines[i - start_idx]
                
            # 写入更新后的mol2文件
            with open(output_file_path, 'w') as f_final_mol2:
                f_final_mol2.writelines(mol2_lines)
            
            print(f"Processed {mol2_file} successfully.")

            # 读取并打印mol2文件的内容
            with open('RESP/'+mol2_file, 'r') as mol2_print:
                mol2_content = mol2_print.read()
                print(mol2_content)

def calculate_capcharge(resp_folder):

    for mol2_file in os.listdir(resp_folder):
        if mol2_file.endswith('.mol2'):
            # 提取文件名前三个字符作为res对象
            res = mol2_file[:3]
            
            # 读取mol2文件内容
            with open(os.path.join(resp_folder, mol2_file), 'r+') as f_mol2:
                mol2 = f_mol2.readlines()
                
                # 定位ATOM所在行
                end = mol2.index('@<TRIPOS>BOND\n') 
                start = mol2.index('@<TRIPOS>ATOM\n')

                ace_cap, nme_cap = 0, 0

                # 计算ACE封端电荷
                for i in range(1 + start, 7 + start):
                    lis = list(filter(None, mol2[i].replace('\n', '').split(' ')))
                    ace_cap += eval(lis[-1])
                
                # 计算NME封端电荷
                for j in range(7 + start, 13 + start):
                    lis = list(filter(None, mol2[j].replace('\n', '').split(' ')))
                    nme_cap += eval(lis[-1])

                # 四舍五入电荷值，以保留六位小数
                ace_cap, nme_cap = round(ace_cap, 6), round(nme_cap, 6)
                
                # 写入封端电荷信息
                charge_info = open(f'{res}_cap.charge', 'w')
                charge_info.write(f'ace_cap: {ace_cap}\n')
                print(f'Sum charge of ACE: {ace_cap}')
                charge_info.write(f'nme_cap: {nme_cap}\n\n')
                print(f'Sum charge of NME: {nme_cap}')
                charge_info.close()
                
                # 删除生成的电荷信息文件
                os.remove(f'{res}_cap.charge')
                
                # N端与ACE封端电荷相加
                N_ncaa_line = mol2[start + 13]
                N_ncaa = list(filter(None, N_ncaa_line.replace('\n', '').split(' ')))
                if N_ncaa[1] == 'N':
                    print(f'Detect N-termini, with original charge {N_ncaa[-1]}')
                    new_N_charge = round(eval(N_ncaa_line[-10:-1]) + ace_cap, 6)
                    new_N_charge_str = f'{new_N_charge:9.6f}'  # 确保电荷值的格式化
                    mol2[start + 13] = N_ncaa_line.replace(N_ncaa_line[-10:-1], str(new_N_charge), 1)
                    print(f'Update N-termini charge with new value {new_N_charge}')
                else:
                    raise Exception("Please check your PDB input, ensure N-CA-C-O order.")

                # C端与NME封端电荷相加
                C_ncaa_line = mol2[start + 15]
                C_ncaa = list(filter(None, C_ncaa_line.replace('\n', '').split(' ')))
                if C_ncaa[1] == 'C':
                    print(f'Detect C-termini, with original charge {C_ncaa[-1]}')
                    new_C_charge = round(eval(C_ncaa_line[-10:-1]) + nme_cap, 6)
                    new_C_charge_str = f'{new_C_charge:9.6f}'  # 确保电荷值的格式化
                    mol2[start + 15] = C_ncaa_line.replace(C_ncaa_line[-9:-1], str(new_C_charge), 1)
                    print(f'Update C-termini charge with new value {new_C_charge}')
                else:
                    raise Exception("Please check your PDB input, ensure N-CA-C-O order.")

def process_params_file(n_value):
    filename = f'{n_value}.params'
    
    # 读取文件内容
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # 查找并处理以“ICOOR_INTERNAL    O ”开头的行
    for i, line in enumerate(lines):
        if line.startswith('ICOOR_INTERNAL    O '):
            parts = line.split()
            if len(parts) == 8:
                dihedral_angle = float(parts[2])
                bond_angle = float(parts[3])
                
                # 检查二面角是否在180°左右，键角是否在60°左右，若不在则对该行进行修正
                if not (179 <= dihedral_angle <= 181 or -181 <= dihedral_angle <= -179) or not (59 <= bond_angle <= 61):
                    new_line = 'ICOOR_INTERNAL    O    179.999969   59.199989    1.231015   C     CA  UPPER\n'
                    lines[i] = new_line
                    break
    
    # 将修改后的内容写回文件
    with open(filename, 'w') as file:
        file.writelines(lines)

    # 读取并打印params文件的内容
    with open(filename, 'r') as params_print:
        params_content = params_print.read()
        print(params_content)

def read_mol2_file(mol2_filepath):
    atom_lines = []
    with open(mol2_filepath, 'r') as f:
        lines = f.readlines()
    
    atom_started = False
    for line in lines:
        if line.startswith('@<TRIPOS>ATOM'):
            atom_started = True
            continue
        if line.startswith('@<TRIPOS>BOND'):
            atom_started = False
            continue
        if atom_started and line.strip():  # only capture non-empty lines between ATOM and BOND sections
            atom_lines.append(line.strip())
    
    return atom_lines[12:]  # Skip the first 12 lines

def read_params_file(params_filepath):
    atom_lines = []
    with open(params_filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_lines.append(line.strip())
    
    return atom_lines

def adjust_charges_to_integer(charges_list,system_charge):
    #计算当前净电荷距离整数电荷的差值
    rounded_charges = [round(charge, 2) for charge in charges_list]

    charge_diff = system_charge - sum(rounded_charges)
    
    if charge_diff == 0:
        return rounded_charges

    # 按照电荷绝对值从大到小排序
    sorted_indices = sorted(range(len(charges_list)), key=lambda i: abs(charges_list[i]), reverse=True)

    # 按电荷从大到小依次进行±0.01的调整，从而在对体系电荷产生最小影响的情况下使电荷值为整数
    for _ in range(abs(int(charge_diff * 100))):  # 需要的调整次数
        for i in sorted_indices:
            if charge_diff > 0:
                rounded_charges[i] += 0.01
                charge_diff -= 0.01
            elif charge_diff < 0:
                rounded_charges[i] -= 0.01
                charge_diff += 0.01
            if round(charge_diff, 2) == 0:
                break
    
    return rounded_charges

def modify_params_file(params_filepath, rounded_charges):
    with open(params_filepath, 'r') as f:
        lines = f.readlines()

    atom_lines = []
    atom_indices = []

    #读取params文件的ATOM行
    for idx, line in enumerate(lines):
        if line.startswith('ATOM'):
            atom_lines.append(line.strip())
            atom_indices.append(idx)

    #检查原子数和电荷列表数是否匹配
    if len(atom_lines) != len(rounded_charges):
        print(f"ATOM lines count: {len(atom_lines)}")
        print(f"Rounded charges count: {len(rounded_charges)}")
        raise ValueError('params 文件中的 ATOM 行数与电荷列表长度不匹配。')

    for i in range(len(atom_lines)):
        original_line = atom_lines[i]
        # 提取新的电荷值
        charge_part = original_line[-5:]  # 提取最后的电荷部分
        new_charge = f"{rounded_charges[i]:.2f}"
        
        # 检查正数电荷值并在前面加空格以保证格式统一
        if float(new_charge) > 0:
            new_charge = f" {new_charge}"
        
        new_line = original_line[:-5] + new_charge + '\n'  # 应用新的电荷值组装新的行内容
        lines[atom_indices[i]] = new_line  # 更新原始文件中的ATOM行

    # 将修改后的 ATOM 行写回文件
    with open(params_filepath, 'w') as f:
        f.writelines(lines)

def update_params_file_with_temps(params_filepath, temps_filepath):
    # 读取params文件中的所有行
    with open(params_filepath, 'r') as f:
        params_lines = f.readlines()

    # 读取temps文件中的所有以ATOM开头的行
    temps_atoms = []
    with open(temps_filepath, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                temps_atoms.append(line.strip())

    # 遍历params文件中的所有以ATOM开头的行，进行替换
    updated_params_lines = []
    for line in params_lines:
        if line.startswith('ATOM'):
            atom_id = line[:18]
            for temps_line in temps_atoms:
                if temps_line[:18] == atom_id:
                    line = temps_line + '\n'
                    break
        updated_params_lines.append(line)

    # 将更新后的行写回params文件
    with open(params_filepath, 'w') as f:
        f.writelines(updated_params_lines)

    # 读取并打印params文件的内容
    with open(params_filepath, 'r') as params_new_print:
        params_new_content = params_new_print.read()
        print(params_new_content)
      
def process_mol2_file(file_path):
    # Read the mol2 file into a list of lines
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Process each line according to the specified conditions
    for i in range(len(lines)):
        line = lines[i]
        if len(line) >= 52 and line[50:52] == "DU" or line[50:52] == "N3":
            if "N" in line[8:10]:
                lines[i] = line[:50] + "N " + line[52:]

    # Write the modified lines back to the original file
    with open(file_path, 'w') as f:
        f.writelines(lines)
        
def process_mol2_with_args(mol2_value):
    # Construct the file name based on the argument
    file_name = f"{mol2_value}_1.mol2"

    # Get the full file path
    folder_path = 'RESP'  # Assuming RESP folder is in the current working directory
    file_path = os.path.join(folder_path, file_name)

    # Process the mol2 file
    process_mol2_file(file_path)
    
    # 读取并打印mol2文件的内容
    with open(file_path, 'r') as mol2_new_print:
        mol2_new_content = mol2_new_print.read()
        print(mol2_new_content)
        
def generate_top(resp_folder):

    # 遍历RESP文件夹中的mol2文件
    for mol2_file in os.listdir(resp_folder):
        if mol2_file.endswith('.mol2'):
            # 提取文件名前三个字符作为res对象
            res = mol2_file[:3]
            global res_rtp
            res_rtp=res

            # 创建存放结果的文件夹
            output_folder = f'./{res}_gromacs_prm'
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            # 执行parmchk2命令生成.mod文件
            mod_file = f'{res}.mod'
            os.system(f'parmchk2 -i {os.path.join(resp_folder, mol2_file)} -f mol2 -o {mod_file}')

            # 生成leapin文件
            leapin_filename = f'{res}_leap.in'
            with open(leapin_filename, 'w+') as leapin:
                leapin.write(f'source /home/rotations/.conda/envs/rdkit/dat/leap/cmd/leaprc.protein.ff19SB\n')
                leapin.write(f'loadamberparams {mod_file}\n')
                leapin.write(f'mol=loadmol2 {os.path.join(resp_folder, mol2_file)}\n')
                leapin.write(f'check mol\n')
                leapin.write(f'saveamberparm mol {res}.prm {res}.crd\n')
                leapin.write('quit\n')

            # 调用tleap
            os.system(f'tleap -f {leapin_filename}')

            # 移动生成的文件到指定文件夹
            generated_files = [f'{res}.prm', f'{res}.crd', mod_file, leapin_filename, 'leap.log']
            for file in generated_files:
                if os.path.exists(file):
                    shutil.move(file, os.path.join(output_folder, file))

            # 切换到输出文件夹进行ACpype操作
            os.chdir(output_folder)
            os.system(f'acpype -p {res}.prm -x {res}.crd -c user -o gmx -a amber')

            # 返回到原始工作目录
            os.chdir('..')

            # 移动生成的GROMACS文件到指定文件夹
            gromacs_files = ['MOL_GMX.gro', 'MOL_GMX.top']
            for file in gromacs_files:
                if os.path.exists(file):
                    shutil.move(file, os.path.join(output_folder, f'{res}.{file.split("_")[-1]}'))

            # 清理中间文件
            os.system(f'rm qout QOUT punch md.mdp esout em.mdp')

    print(f'\nFZ-wang reminds you: The GROMACS top files have been generated in the folder "gromacs_prm"!\n')
    

def generate_rtp(resp_file):
    res = res_rtp  # 残基的名称或标识
    ignore_atoms = ["ACE", "NME", "linker"]  # 不需要包含在 .rtp 文件中的原子名称列表
    # 打开并读取 .top 文件
    with open(resp_file, 'r') as f_top:
        top = f_top.readlines()

    # 确定各个部分的起始和结束行索引
    start_atom = top.index('[ atoms ]\n') + 2
    start_bond = top.index('[ bonds ]\n') + 2
    end_atom = top.index('[ bonds ]\n') - 1
    end_bond = top.index('[ pairs ]\n') - 1
    start_angle = top.index('[ angles ]\n') + 2
    end_angle = top.index('[ dihedrals ] ; propers\n') - 1
    start_dihedral = top.index('[ dihedrals ] ; propers\n') + 3
    end_dihedral = top.index('[ dihedrals ] ; impropers\n') - 1
    start_improper = top.index('[ dihedrals ] ; impropers\n') + 3
    end_improper = top.index('[ system ]\n') - 1

    # 初始化 RTP 列表
    rtp_list = []
    rtp_list.append(f'[ {res} ]\n')  # 残基条目
    rtp_list.append(' [ atoms ]\n')  # atoms

    # 定义 include_ffparm 函数，用于排除不需要的行
    def include_ffparm(atom_nums, atom_names):
        for i in atom_nums:
            if int(i) < 13:
                return False
        for j in atom_names:
            if j in ignore_atoms:
                return False
        return True
    
    # 处理 atoms 项
    for i in range(start_atom, end_atom):
        line = top[i]
        atom = list(filter(None, line.replace('\n', '').split()))

        atom_nums, atom_names = [atom[0]], [atom[4]]
        atom_name, atom_type, atom_charge, atom_num = atom[4], atom[1], atom[6], int(atom[0])-12
        
        if include_ffparm(atom_nums, atom_names):
#           print(f'processing atom {atom_name}')
            rtp_list.append(f'    {atom_name:>4}   {atom_type:>2}    {atom_charge:>9}    {atom_num:>2}\n')

    rtp_list.append('\n [ bonds ]\n')

    # 处理 bonds 项
    for j in range(start_bond, end_bond):
        line = top[j]
        bond = list(filter(None, line.replace('\n', '').split()))

        atom_nums, atom_names = bond[0:2], [bond[-3], bond[-1]]
        if include_ffparm(atom_nums, atom_names):
#           print(f'processing bond {bond[-3]}-{bond[-1]}')
            rtp_list.append(f'    {bond[-3]:>4}   {bond[-1]:<4}  {bond[3]}    {bond[4]}\n')

    rtp_list.append(f'    {"-C":>4}   {"N":<4}  1.3790e-01    3.5782e+05\n')
    rtp_list.append('\n [ angles ]\n')

    # 处理 angles 项
    for k in range(start_angle, end_angle):
        line = top[k]
        angle = list(filter(None, line.replace('\n', '').split()))

        atom_nums, atom_names = angle[0:3], [angle[-5], angle[-3], angle[-1]]
        if include_ffparm(atom_nums, atom_names):
#           print(f'processing angle {angle[-5]}-{angle[-3]}-{angle[-1]}')
            rtp_list.append(f'    {angle[-5]:>4}   {angle[-3]:>4}    {angle[-1]:<4}  {angle[4]}   {angle[5]}\n')

    rtp_list.append('\n [ dihedrals ] ; propers\n')

    # 处理 dihedrals 项
    for l in range(start_dihedral, end_dihedral):
        line = top[l].replace('-', ' ')
        dihedral = list(filter(None, line.replace('\n', '').split()))

        atom_nums, atom_names = dihedral[0:4], [dihedral[-4], dihedral[-3], dihedral[-2], dihedral[-1]]
        if include_ffparm(atom_nums, atom_names):
#           print(f'processing dihedrals proper {dihedral[-4]}-{dihedral[-3]}-{dihedral[-2]}-{dihedral[-1]}')
            rtp_list.append(f'    {dihedral[-4]:>4}   {dihedral[-3]:>4}   {dihedral[-2]:>4}   {dihedral[-1]:<4}  {dihedral[5]:>6}   {dihedral[6]:>8}   {dihedral[7]}\n')

    rtp_list.append('\n [ dihedrals ] ; impropers\n')

    # 处理 impropers 项
    for m in range(start_improper, end_improper):
        line = top[m].replace('-', ' ')
        improper = list(filter(None, line.replace('\n', '').split()))

        atom_nums, atom_names = improper[0:4], [improper[-4], improper[-3], improper[-2], improper[-1]]
        if include_ffparm(atom_nums, atom_names):
#           print(f'processing dihedrals impropers {improper[-4]}-{improper[-3]}-{improper[-2]}-{improper[-1]}')
            rtp_list.append(f'    {improper[-4]:>4}   {improper[-3]:>4}   {improper[-2]:>4}   {improper[-1]:<4}  {improper[5]:>6}   {improper[6]:>8}   {improper[7]}\n')

    rtp_list.append('    -C    CA     N     H  180.00   4.60240   2\n    CA    +N     C     O  180.00   4.60240   2\n')

    # 写入生成的 .rtp 文件
    with open(f'{res_rtp}.rtp', 'w') as f_rtp:
        for line in rtp_list:
            f_rtp.write(line)

    # 读取并打印rtp文件的内容
    with open(f'{res_rtp}.rtp', 'r') as rtp_print:
        rtp_content = rtp_print.read()
        print(rtp_content)
        
def main(input_file, names_list, clean):
    # 确保所有需要的文件夹存在或根据需要创建它们
    required_dirs = ['pdb_files', 'PDB_rearranged', 'GJF', 'mol']
 
    
    gjf_folder = 'GJF'
    resp_folder = 'RESP'
    
    for dir_name in required_dirs:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            print(f"Created directory: {dir_name}")
 
    # 读取输入的SMILES文件
    
    # 计算体系电荷
    system_charge = calculate_system_charge(input_file)
    print(f"System charge: {system_charge}")
    
    with open(input_file, 'r') as f:
        smiles_list = [line.strip() for line in f.readlines()]
 
    if len(names_list) != len(smiles_list):
        raise ValueError("The number of provided names does not match the number of SMILES strings.")
    
    # 处理SMILES
    processed_smiles_list = []
    for smiles in smiles_list:
        try:
            processed_smiles = process_smiles(smiles)
            processed_smiles_list.append(processed_smiles)
        except ValueError as e:
            print(e)
            continue
 
    # 将处理后的SMILES转换为PDB
    output_dir = 'pdb_files'
    generated_files = smiles_to_pdb(processed_smiles_list, names_list, output_dir)
 
    # 重新排序并处理PDB文件
    output_rearranged_dir = 'PDB_rearranged'
    process_files(generated_files, names_list, output_rearranged_dir)
 
    # 调整二面角
    adjusted_pdb_files = [os.path.join(output_rearranged_dir, f'{name}_{i+1}.pdb') for i, name in enumerate(names_list)]
    adjust_dihedrals_in_pdb_files(adjusted_pdb_files)
 
    # 删除PDB_rearranged文件夹中所有PDB文件的CONECT行
    for pdb_file in os.listdir(output_rearranged_dir):
        if pdb_file.endswith(".pdb"):
            pdb_filepath = os.path.join(output_rearranged_dir, pdb_file)
            remove_conect_lines_from_pdb(pdb_filepath)
 
    # 生成Gaussian输入文件
    gaussian_input_dir = 'GJF'
    generate_gaussian_input(adjusted_pdb_files, gaussian_input_dir,system_charge)
    
    # 等待所有Gaussian作业完成
    all_jobs_completed = False
    while not all_jobs_completed:
        all_jobs_completed = True
        for log_file in os.listdir(gaussian_input_dir):
            if log_file.endswith('.log'):
                log_filepath = os.path.join(gaussian_input_dir, log_file)
                if not is_gaussian_job_completed(log_filepath):
                    all_jobs_completed = False
                    break
        if not all_jobs_completed:
            time.sleep(60)  # 等待1分钟后再次检查未完成的作业
    
    # 所有Gaussian作业完成后，处理log文件并生成mol文件
    process_log_files(input_file, gaussian_input_dir)
    
        # 获取生成的mol文件列表
    mol_files = [f for f in os.listdir('mol') if f.endswith('.mol')]
    
    # 遍历每个mol文件，转化为params文件并存放到对应的文件夹中
    for mol_file in mol_files:
        base_name = os.path.splitext(mol_file)[0]
        name = base_name.split('_')[0]
        mol_filepath = os.path.join('mol', mol_file)
        molfile_to_params(mol_filepath, name)
        molfile_to_params_temps(mol_filepath, name)
    
    # 确保RESP文件夹存在
    if not os.path.exists(resp_folder):
        os.makedirs(resp_folder)
    
    # 生成resp文件
    for filename in os.listdir(gjf_folder):
        if filename.endswith('.log'):
            res = os.path.splitext(filename)[0]  # Get the base name without extension
            resp = os.path.join(resp_folder, res)  # RESP file path
    
            # Copy log file to current working directory
            shutil.copyfile(os.path.join(gjf_folder, filename), f'{res}.log')
    
            # Run the generate_opt function
            generate_opt(res, resp)
    
            # Move the log file back to the RESP folder
            shutil.move(f'{res}.log', os.path.join(resp_folder, f'{res}.log'))
    
    # 处理原子类型调整
    atom_type_adjust(resp_folder)
    
    # 计算capcharge
    calculate_capcharge(resp_folder)
    
    n_value = args.names[0]
    # 处理params文件
    process_params_file(n_value)
    
    mol2_filepath = os.path.join('RESP', f'{n_value}_1.mol2')
    params_filepath = os.path.join(f'{n_value}_temps.params')
 
    # Step 1: Calculate system charge
    system_charge = calculate_system_charge('input.smiles')
    print(f"System charge: {system_charge}")
 
    # Step 2: Read and process mol2 file
    atom_lines = read_mol2_file(mol2_filepath)
    charges_list = []
    for line in atom_lines:
        parts = line.split()
        charge = float(parts[-1])
        charges_list.append(charge)
    
    # Step 3: Adjust charges to ensure they sum to an integer
    rounded_charges = adjust_charges_to_integer(charges_list,system_charge)
    print(f"Charges list after adjustment: {rounded_charges}")
 
    # Step 4: Modify params file
    modify_params_file(params_filepath, rounded_charges)
    print(f"Modified {params_filepath} with rounded charges.")
 
    # Step 5: Update params file with temps params
    temps_filepath = os.path.join(f'{n_value}.params')
    update_params_file_with_temps(temps_filepath, params_filepath)
    print(f"Updated {params_filepath} with data from {temps_filepath}")
    
    process_mol2_with_args(n_value)
    
    # 生成top文件
    generate_top(resp_folder)
    
    resp = f"{res_rtp}_gromacs_prm/MOL.amb2gmx/MOL_GMX.top"
 
    # 生成rtp文件
    generate_rtp(resp)
    
    
    # 删除中间文件夹
    if clean == 1:
        shutil.rmtree('pdb_files', ignore_errors=True)
        shutil.rmtree('PDB_rearranged', ignore_errors=True)
        shutil.rmtree('GJF', ignore_errors=True)
        shutil.rmtree('mol', ignore_errors=True)
        shutil.rmtree('RESP', ignore_errors=True)
        os.remove('molecule.chk')
        os.remove('fort.7')
        print("Intermediate folders deleted.")
        

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process SMILES and convert to PDB files.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input SMILES file.')
    parser.add_argument('-n', '--names', required=True, nargs='+', help='Names for each SMILES (UAA if not specified).')
    parser.add_argument('-c', '--clean', type=int, choices=[0, 1], default=1, help='Clean intermediate folders: 1 to clean (default), 0 to keep.')
    args = parser.parse_args()
    main(args.input, args.names, args.clean)
