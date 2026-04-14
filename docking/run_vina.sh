#!/bin/bash
#SBATCH --job-name=petase_vina
#SBATCH --partition=cpudebug
#SBATCH --qos=cpudebug
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/work/bio/shuleihe23/docking/logs/%x_%j.out
#SBATCH --error=/gpfs/work/bio/shuleihe23/docking/logs/%x_%j.err

set -euo pipefail

echo "[INFO] Job started at: $(date)"
echo "[INFO] Running on node: $(hostname)"

module purge
module load autodock-vina/1.2.5-gcc-8.5.0-ugdue5u

WORKDIR="/gpfs/work/bio/shuleihe23/docking"
mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

echo "[INFO] Working directory: $(pwd)"
echo "[INFO] Using vina at: $(which vina)"

vina \
  --config config.txt \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --out out.pdbqt \

echo "[INFO] Job finished at: $(date)"