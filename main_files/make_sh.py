text ="""#!/bin/bash
#SBATCH --job-name={}
#SBATCH --partition=q24,q20,q16
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=machri@phys.au.dk

echo "========= Job started  at `date` =========="

source /home/machri/envs/python3.6.3/bin/activate

cd $SLURM_SUBMIT_DIR
cp *.py /scratch/$SLURM_JOB_ID

python main.py
cp *.out $SLURM_SUBMIT_DIR 
cp *.xyz $SLURM_SUBMIT_DIR 

echo "========= Job finished at `date` =========="
"""


main_name = 'test'

num_reps = 1
num_files = 1

jobs = []
for jj in range(0, num_files):
    job_name = main_name+'_{}_{}'.format(num_reps, jj)
    jobs.append(job_name)
    with open(job_name+'.sh', 'w') as f:
        print(text.format(job_name, num_reps, jj), file=f)


with open('run_file.sh', 'w') as f:
    print('#!/bin/bash', file=f)

    for job_name in jobs:
        print('sbatch '+job_name+'.sh', file=f)




