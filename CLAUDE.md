# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

gen2ntuple creates ROOT ntuples from MicroBooNE LAr TPC reconstruction data. It processes larflowreco and merged_dlreco files to produce flat data structures for physics analysis, incorporating deep learning particle identification and Monte Carlo truth information.

## Essential Commands

### Initial Setup (One-time)
```bash
# Compile the C++ fiducial volume library
cd helpers/
source compile_wirecell_fiducial_volume.sh
cd ..
```

### Running the Ntuple Maker
```bash
# Basic usage
python make_dlgen2_flat_ntuples.py \
  -f <larflowreco_files> \
  -t <merged_dlreco_files> \
  -w <weight_file> \
  -m <larpid_model_path> \
  -o <output_file> \
  [-mc] [-ana] [-nkp]

# Example for MC overlay data
python make_dlgen2_flat_ntuples.py \
  -f larflow_file.root \
  -t merged_dlreco.txt \
  -w event_weighting/weights_forCV_v48_Sep24_bnb_nu_run1.pkl \
  -m /path/to/ResNet34_recoProng_5class_epoch20.pt \
  -o output_ntuple.root \
  -mc
```

### Testing/Verification
```bash
# Run example analysis
python example_ntuple_analysis_script.py -i ntuple.root -o histograms.root

# Verify photon energy deposition (specific test)
python tests/verify_photon_edep.py
```

### Batch Processing on Tufts Cluster
```bash
# Submit SLURM job
sbatch tufts_submit_ntuple_job_example.sh
```

## High-Level Architecture

### Core Workflow
1. **Input Stage**: Reads larflowreco (reconstruction) and merged_dlreco (truth) files
2. **Processing Stage**: 
   - Identifies neutrino vertices from keypoint data
   - Reconstructs tracks and showers
   - Runs LArPID CNN for particle classification
   - Calculates physics variables (energy, containment, angles)
3. **Output Stage**: Writes flat ROOT ntuples with EventTree and potTree

### Key Components

- **make_dlgen2_flat_ntuples.py**: Main entry point, orchestrates the workflow
- **helpers/larflowreco_ana_funcs.py**: Core analysis functions for reconstruction
- **helpers/select_nu_vertex.py**: Neutrino vertex selection algorithms
- **helpers/intime_vertex.py**: In-time cosmic rejection
- **helpers/pionEnergyEstimator.py**: Pion kinetic energy reconstruction
- **event_weighting/**: Pre-computed cross-section weights for different runs

### Data Flow
```
larflowreco files → Vertex Selection → Track/Shower Reconstruction → 
LArPID Classification → Energy Calculation → ROOT Ntuple Output
```

### Key Dependencies
- ROOT (via PyROOT)
- PyTorch (for LArPID network)
- ubdl environment (larlite, larcv, larflow, ublarcvapp)
- Compiled C++ fiducial volume library

### Important Notes
- The project uses Wire Cell fiducial volume definition (>3cm from SCE-corrected edges)
- Event weights are crucial for proper POT scaling in MC samples
- LArPID model path must point to a trained ResNet34 5-class model
- Output ntuples contain 99+ physics variables documented in README.md