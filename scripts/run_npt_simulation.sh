#!/usr/bin/env bash

#SBATCH --job-name=melting_quenching_Au_doped
#SBATCH --partition=cenvalarc.gpu
#SBATCH --constraint=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=2-23:59:00

#SBATCH --output=job.o%j
#SBATCH --error=job.e%j

###PreProcessing###
module purge
module load anaconda3
module load fftw
module load cuda
# Try to load a specific older OpenMPI version that's more stable
module load openmpi/4.1.1-gcc-8.4.1

# Initialize conda in the current shell
eval "$(conda shell.bash hook)"
# Activate conda environment
echo "Activating conda environment: fair"
conda activate fair

echo "========================================"
echo "MoS2 NPT Simulation Runner"
echo "========================================"

# Check if we're in the right directory
if [ ! -f "src/mos2_npt_simulation.py" ]; then
    echo "Error: Please run this script from the Doped_MoS2_density directory"
    exit 1
fi

# Create results directory if it doesn't exist
mkdir -p results
cd results

# Copy the simulation script to results directory for this run
cp ../src/mos2_npt_simulation.py ../src/surfaces.py .

# Run the simulation
echo "Starting MoS2 NPT simulation..."
echo "Simulation started at: $(date)"

# Run with error handling
python mos2_npt_simulation.py 2>&1 | tee simulation_output.log

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully at: $(date)"
    echo "Results saved in: $(pwd)"
else
    echo "Simulation failed. Check simulation_output.log for details."
    exit 1
fi

echo "========================================"
echo "Simulation Complete"
echo "========================================"
