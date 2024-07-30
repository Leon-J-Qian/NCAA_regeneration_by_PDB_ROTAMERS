#!/bin/bash
#SBATCH -J Mak012553
#SBATCH -p cn-long
#SBATCH -N 1 
#SBATCH -o Mak012553_%j.out
#SBATCH -e Mak012553_%j.err
#SBATCH --no-requeue
#SBATCH -A chuwang_g1
#SBATCH --qos=chuwangcnl
#SBATCH --ntasks-per-node=20
hosts=`scontrol show hostname $SLURM_JOB_NODELIST` ;echo $hosts
. /apps/source/intel-2018.sh
cat ./*.list | parallel -j 20