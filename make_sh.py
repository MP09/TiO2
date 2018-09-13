import sys
import os
import glob

text ="""#!/bin/bash
#SBATCH --job-name={}
#SBATCH --partition=q24,q20,q16
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=out/{}

echo "========= Job started  at `date` =========="
echo $SLURM_SUBMIT_DIR


source /home/machri/envs/python3.6.3/bin/activate
cd $SLURM_SUBMIT_DIR
cp *.py /scratch/$SLURM_JOB_ID

cd /scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=${}

python bh_runner.py {} {} {} {} {}
cp *.traj $SLURM_SUBMIT_DIR/trajs/
cp *.npy  $SLURM_SUBMIT_DIR/npys/
cp *.out  $SLURM_SUBMIT_DIR/out/
cp *.fail $SLURM_SUBMIT_DIR/fail/
echo "========= Job finished at `date` =========="
"""

main_name = ''
spc = r'''{SLURM_CPUS_PER_TASK:-1}'''

# Inputs:
ML_setting = int(sys.argv[1])
pertub_setting = int(sys.argv[2])
num_files = int(sys.argv[3])
ncalc = int(sys.argv[4])
max_iter = int(sys.argv[5])

# Pertubation settings: 
if pertub_setting == 0:
    folder1 = 'FW/'
    extra_name = '_FW'
    label = 'fireworks'
elif pertub_setting == 1:
    folder1 = 'RC/'
    extra_name = '_RC'
    label = 'random_center'
elif pertub_setting == 2:
    folder1 = 'RP/'
    extra_name = '_RP'
    label = 'random_position'

# ML settings:
if ML_setting == 1:
    ml_name = 'ML'
    folder1 = 'ML_'+folder1  
else:
    ml_name = ''


folder2 = 'bash/'
folder = folder1+folder2
# Make folders:
if not os.path.exists(folder):
    os.makedirs(folder)
    os.mkdir(folder1+'trajs')
    os.mkdir(folder1+'out')
    os.mkdir(folder1+'fail')
    os.mkdir(folder1+'npys')

# Write bash files:
jobs = []
for jj in range(0, num_files):
    job_name = ml_name+main_name+extra_name+'_{}'.format(jj)
    out_name = ml_name[0:-1]+main_name+extra_name+'_{}.out'.format(jj)
    jobs.append(job_name)
    with open(folder+job_name+'.sh', 'w') as f:
        ninit = ncalc*jj
        print(text.format(job_name, out_name, spc, ncalc, ninit, pertub_setting, ML_setting, max_iter), file=f)

# Write run_file
run_file = ml_name+'run_file.sh'

run_file = folder1 + run_file
with open(run_file, 'w') as f:
    print('#!/bin/bash', file=f)
    for job_name in jobs:
        print('sbatch '+folder2+job_name+'.sh', file=f)


# Information out:
from pertubations import pertubation_dict

ml_dict = {0:'without', 1:'with'}
print('='*50)
pertubation = pertubation_dict[pertub_setting][0]

print('Made {} files with the {} pertubation {} autobagging'.format(num_files, pertubation, ml_dict[ML_setting]))
print('='*50)




