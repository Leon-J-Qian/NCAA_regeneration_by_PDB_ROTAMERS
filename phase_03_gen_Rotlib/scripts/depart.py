import sys

ncaa = sys.argv[1]
num_chi = int(sys.argv[2])
rotwell_value = sys.argv[3:]

rosetta_path = '~/rosetta/rosetta_src_2015.19.57819_bundle/main/source/bin/MakeRotLib.mpi.linuxiccrelease'

def generate_chi_wells(num_chi, rotwell_value):
    assert len(rotwell_value) == num_chi

    rotwells_table = {'2':[' 0', '180'], '3':[' 60', '180', '300'],'S':[' 90', '180', '270'], '6':['  0', ' 60', '120', '180', '240', '300'], '12':['  0', ' 30', ' 60', ' 90', '120', '150', '180', '210', '240', '270', '300', '330']}
    rotwells = []
    for i in rotwell_value:
        rotwells.append(rotwells_table[i])

    return rotwells

def dfs(arr, i, j, cur_seq, result):
    cur_seq.append((arr[i][j], j+1))
    if i!=len(arr)-1:
         dfs(arr,i+1,0,cur_seq, result)
    else:
         result.append(cur_seq)

    if j!=len(arr[i])-1:
        new_seq = cur_seq[:i]
        dfs(arr,i,j+1,new_seq, result)
    
    return result

f_job = open(f'MRL_job_{ncaa}.list', 'w')

for phi in range(-180,190,10):
    for psi in range(-180,190,10):
        f = open(f'{ncaa}_{phi}_{psi}.in','w')
        f.write(f'AA_NAME {ncaa}\n')
        f.write(f'OMG_RANGE 180 180 1\nPHI_RANGE {phi} {phi} 1\nPSI_RANGE {psi} {psi} 1\nEPS_RANGE 180 180 1\n')
        f.write(f'NUM_CHI {num_chi}\nNUM_BB 2\n')
        for i in range(num_chi):
            f.write(f'CHI_RANGE {i+1} 0  330  30\n')
        rotwells = generate_chi_wells(num_chi, rotwell_value)
        results = dfs(rotwells, 0, 0, [], [])
        for i in results:
            cen = 'CENTROID '
            for j in i:
                cen += f'{j[0]} {j[1]:>2} '
            f.write(cen+'\n')
        f.close()
        f_job.write(f'{rosetta_path} -options_file {ncaa}_{phi}_{psi}.in -extra_res_fa {ncaa}.params -mute all -score:weights mm_std ; rm -rf {ncaa}*mrllog \n')
