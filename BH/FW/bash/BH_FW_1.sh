#!/bin/bash
#SBATCH --job-name=BH_FW_1
#SBATCH --partition=q24,q20,q16
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --output=out/BH_1.out

echo "========= Job started  at `date` =========="
echo $SLURM_SUBMIT_DIR


source /home/machri/envs/python3.6.3/bin/activate
cd $SLURM_SUBMIT_DIR
cp *.py /scratch/$SLURM_JOB_ID

cd /scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python bh_run.py 1 2000 0 fireworks
cp *.traj $SLURM_SUBMIT_DIR/trajs/
cp *.npy  $SLURM_SUBMIT_DIR/npys/
cp *.out  $SLURM_SUBMIT_DIR/out/
cp *.fail $SLURM_SUBMIT_DIR/fail/
echo "========= Job finished at `date` =========="

