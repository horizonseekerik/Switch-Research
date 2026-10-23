#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AUTOMATED MINIFORGE & PYMEEP INSTALLER FOR AZURE VM
# (Requires ZERO root/sudo privileges - installs cleanly into $HOME/miniforge3)
# ==============================================================================
set -e

echo "======================================================================"
echo "  PROJECT JANUS: PYMEEP FDTD INSTALLER (AZURE CLOUD / VM)"
echo "======================================================================"

# 1. Download and install Miniforge if not already installed
if [ ! -d "$HOME/miniforge3" ]; then
    echo "[1/3] Downloading Miniforge installer..."
    curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    
    echo "[2/3] Installing Miniforge into $HOME/miniforge3..."
    bash Miniforge3.sh -b -p "$HOME/miniforge3"
    rm -f Miniforge3.sh
    
    "$HOME/miniforge3/bin/conda" init bash
else
    echo "[INFO] Miniforge is already installed at $HOME/miniforge3"
fi

# 2. Activate Conda environment in current subshell
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate base

# 3. Install PyMEEP from conda-forge
echo "[3/3] Installing pymeep and dependencies from conda-forge (this may take 2-3 mins)..."
conda install -c conda-forge pymeep python=3.11 -y

# 4. Verify installation
echo ""
echo "======================================================================"
echo "  VERIFYING PYMEEP FDTD INSTALLATION"
echo "======================================================================"
python -c "import meep as mp; print(f'SUCCESS! PyMEEP version {mp.__version__} is active and ready for 3D Maxwell simulations.')"

echo ""
echo "======================================================================"
echo "  Ready to run 3D Full-Wave Maxwell FDTD simulations!"
echo "  Command: python sub2um_switch_3d_meep.py --state both --resolution 25"
echo "======================================================================"
