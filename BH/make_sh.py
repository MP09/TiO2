import sys
import os
import glob

text ="""#!/bin/bash
#SBATCH --job-name={}
#SBATCH --partition=q24,q20,q16
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --output=out/{}

echo "========= Job started  at `date` =========="
echo $SLURM_SUBMIT_DIR


source /home/machri/envs/python3.6.3/bin/activate
cd $SLURM_SUBMIT_DIR
cp *.py /scratch/$SLURM_JOB_ID

cd /scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=${}

python bh_run.py {} {} {} {}
cp *.traj $SLURM_SUBMIT_DIR/trajs/
cp *.npy  $SLURM_SUBMIT_DIR/npys/
cp *.out  $SLURM_SUBMIT_DIR/out/
cp *.fail $SLURM_SUBMIT_DIR/fail/
echo "========= Job finished at `date` =========="
"""

main_name = 'BH'
spc = r'''{SLURM_CPUS_PER_TASK:-1}'''
num_iter = 2000

ML_setting = int(sys.argv[1])
pertub_setting = int(sys.argv[3])

# Some settings: 
if pertub_setting == 0:
    folder1 = 'FW/'
    extra_name = '_FW'
    label = 'fireworks'
elif pertub_setting == 1:
    folder1 = 'RC/'
    extra_name = '_RC'
    label = 'random_center'


folder2 = 'bash/'
folder = folder1+folder2

# Make folders:
if not os.path.exists(folder):
    os.makedirs(folder)
    os.mkdir(folder1+'trajs')
    os.mkdir(folder1+'out')
    os.mkdir(folder1+'fail')
    os.mkdir(folder1+'npys')

# Settings passed to python:
num_reps = 1
num_files = int(sys.argv[2])

# Removes old files:
#a = glob.glob(folder+'*.sh')
#for f in a:
#    os.remove(f)

# Write bash files:
jobs = []
for jj in range(0, num_files):
    job_name = main_name+extra_name+'_{}'.format(jj)
    out_name = main_name+'_{}.out'.format(jj)
    if ML_setting == 1:
        job_name = 'ML_'+job_name
        out_name = 'ML_'+out_name

    jobs.append(job_name)
    with open(folder+job_name+'.sh', 'w') as f:
        print(text.format(job_name, out_name, spc, jj, num_iter, ML_setting, label), file=f)

# Write run_file
run_file = 'run_file.sh'
if ML_setting == 1:
    run_file = 'ML_'+run_file

run_file = folder1 + run_file
with open(run_file, 'w') as f:
    print('#!/bin/bash', file=f)
    for job_name in jobs:
        print('sbatch '+folder2+job_name+'.sh', file=f)

# Write combined run file:
if ML_setting == 1:
    other_file = folder1+'run_file.sh'
else:
    other_file = folder1+'ML_run_file.sh'

files = [run_file, other_file]
combined_run_file = folder1+'run.sh'
with open(combined_run_file, 'w') as f:
    print('#!/bin/bash', file=f)
    for file in files:
        if os.path.exists(file):
            with open(file, 'r') as ff:
                lines = ff.readlines()
                for line in lines[1::]:
                    print(line[0:-1], file=f)


# Information out:
pertubation_dict = {0:'"fireworks"', 1:'"random center"'}
ml_dict = {0:'without', 1:'with'}
print('='*50)
pertubation = pertubation_dict[pertub_setting]

print('Made {} files with the {} pertubation {} autobagging'.format(sys.argv[2], pertubation, ml_dict[ML_setting]))
print('='*50)




