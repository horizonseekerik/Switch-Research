#!/bin/bash
# =============================================================================
# PROJECT JANUS: AZURE UBUNTU 22.04 / 24.04 MEEP 3D FDTD DEPLOYMENT SCRIPT
# =============================================================================
# Turnkey setup for running 4-vCPU parallel 3D Meep FDTD simulations on Azure Free Tier
#
# Usage:
#   chmod +x setup_azure_meep.sh
#   ./setup_azure_meep.sh
# =============================================================================

set -e

echo "==============================================================================="
echo "PROVISIONING PARALLEL MEEP 3D FDTD ENVIRONMENT (AZURE 4-vCPU)"
echo "==============================================================================="

# 1. Update OS and system toolchains
sudo apt-get update -y
sudo apt-get install -y build-essential git wget htop libopenmpi-dev openmpi-bin libhdf5-openmpi-dev

# 2. Install Miniconda if not present
if ! command -v conda &> /dev/null; then
    echo "[1/4] Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    conda init bash
else
    echo "[1/4] Conda already installed."
fi

# 3. Create optimized pymeep environment from conda-forge
echo "[2/4] Installing pymeep-parallel with OpenMPI and MPB..."
conda create -n meep_env -c conda-forge pymeep-parallel python=3.10 numpy scipy matplotlib h5py -y

# 4. Activate environment
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate meep_env

echo "[3/4] Validating parallel Meep installation..."
python3 -c "import meep as mp; print('Meep Version:', mp.__version__)"

# 5. Launch Parallel 3D FDTD Simulation
echo "[4/4] Executing 4-vCPU Parallel 3D Switch Verification..."
python3 run_optimized_meep_3d_switch.py --resolution 20 --num-cpus 4 --cycles 0

echo "==============================================================================="
echo "SIMULATION COMPLETED SUCCESSFULLY!"
echo "To run 100M cycle degraded state verification:"
echo "  python3 run_optimized_meep_3d_switch.py --resolution 20 --num-cpus 4 --cycles 100M"
echo "==============================================================================="
