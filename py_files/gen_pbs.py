# generate pbs scripts for katana
# usage
#   python gen_pbs.py m_min m_max
import sys
m_min  = int(sys.argv[1])
m_max  = int(sys.argv[2])
print("#!/bin/bash ")
print("#PBS -N computeE ") 
print("#PBS -l select=1:ncpus=2:mem=24gb ")
print("#PBS -l walltime=200:00:00 ")
print("#PBS -m ae ")
print("#PBS -M qlegia@unsw.edu.au ")
print("#PBS -J %d-%d \n"%(m_min,m_max))
print("module add python/3.6.5 ")
print("cd /home/z9701564/CMB_Pseudo_Cell ")
print("echo “I am now working on job ${PBS_ARRAY_INDEX}” ")
print("python3 comp_matEv2.py 1000 ${PBS_ARRAY_INDEX} 200 ")
