# Flash Prediction

Here we study the addition of the optical flash prediction to the set of variables we can use to select different neutrino interactions.

Eventually, this quantity should be made as part of the Lantern reconstruction.

But to avoid reprocessing that stage, we can use the current info on the neutrino candidates in the KPSRecoManagerTree along with the input dlmerged files
to calculate the flash predictions and study their effectiveness.

## Overview

This project provides a C++ executable that calculates flash predictions for each neutrino candidate within an event. The predictions are saved into a ROOT tree that can be used as a friend tree with existing ntuples.

## Features

- Calculates predicted flash for each neutrino vertex candidate
- Two prediction modes: using all particles or primary particles only
- Computes Sinkhorn divergence between predicted and observed flashes
- Calculates PE differences and ratios
- Outputs ROOT trees compatible with friend tree analysis

## Building

You must be running inside a container.

One option at Tufts:
  ```
  /cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/u20.04_cu111_cudnn8_torch1.9.0_minkowski_npm.sif
  ```

At Fermilab:
  ```
  /cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod/
  ```

1. **Setup environment** (required before building):
   ```bash
   cd /path/to/ubdl
   source setenv_py3.sh
   source configure.sh
   ```

2. **Build the executable**:
   ```bash
   cd studies/flashprediction
   ./build.sh
   ```

   The executable will be created at: `build/flashprediction/calculate_flash_predictions`

## Usage

### Running the executable

```bash
./build/flashprediction/calculate_flash_predictions \
  -d dlmerged_file.root \
  -r reco_file.root \
  -o output_predictions.root \
  [-n num_entries] \
  [-s start_entry] \
  [-t adc_threshold] \
  [-v] \
  [-tb] \
  [-mc]
```

**Parameters:**
- `-d, --dlmerged`: Input dlmerged ROOT file (contains ADC images and opflash)
- `-r, --reco`: Input reco file with KPSRecoManagerTree
- `-o, --output`: Output ROOT file for flash predictions
- `-n, --num-entries`: Number of entries to process (default: all)
- `-s, --start-entry`: Starting entry (default: 0)
- `-t, --threshold`: ADC threshold (default: 10.0)
- `-v, --verbose`: Enable verbose output
- `-tb, --tickbackward`: Use tick backward direction
- `-mc, --mc`: Enable MC mode (calculate distance to true neutrino vertex)

### Submitting batch jobs

We proceed in a few steps.

First, to launch the jobs use `tufts_submit_flash_prediction.sh`.

You'll need to configure the jobs to run for a sample. Example:
```
# mcc9_v28_wctagger_bnb5e19 : run 1 open data
# number of files: 11681
NFILES=40
# number of jobs=293
SAMPLENAME=mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune
RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune.txt
TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v28_wctagger_bnb5e19.txt
MCFLAG=""
```

Note that we refer to filelists that help index the source image files ("merged_dlreco") and result of the lantern reco files ("kpsrecoana").

Next, we use ROOT's hadd utility to merge the files made by each job.
  ```
  python3 run_hadd.py
  ```

Next, we check to see if these files are aligned properly with core ntuple files (if we've retained the "sub"-ntuple files.
  ```
  python3 match_ntuple_files.py
  ```
If all sub-ntuple files align, it'll use hadd to merged these into a final flash prediction file for the sample.

Finally, we use a script to check to see that each event in the ntuple and flash prediction tree
comes from the same event indexed by (run, subrun, event).


### Combining outputs

```bash
python3 combine_flash_predictions.py \
  -i output1.root output2.root ... \
  -o combined_friend_tree.root \
  [-n original_ntuple.root]
```

## Output Structure

The output stores all vertex candidates per event in `std::vector` branches for perfect ntuple alignment:

**Event-level branches:**
- `entry`, `run`, `subrun`, `event`: Event identification
- `n_vertices`: Number of vertex candidates in this event
- `has_vertices`, `has_flash`: Event-level status flags
- `obs_total_pe`, `obs_time`: Observed flash info (same for whole event)
- `obs_pe_per_pmt`: Vector of observed PE per PMT (size 32)

**Vertex-level vector branches (one element per vertex candidate):**
- `pred_total_pe_all`: Vector of total predicted PE using all particles
- `pred_total_pe_primary`: Vector of total predicted PE using primary particles only
- `pred_pe_per_pmt_all`: Vector of vectors (2D) for per-PMT predictions (all particles)
- `pred_pe_per_pmt_primary`: Vector of vectors (2D) for per-PMT predictions (primary only)
- `n_tracks_all`, `n_showers_all`: Vectors of track/shower counts (all particles)
- `n_tracks_primary`, `n_showers_primary`: Vectors of track/shower counts (primary only)
- `sinkhorn_div_all`: Vector of vectors for Sinkhorn divergence (all particles)
- `sinkhorn_div_primary`: Vector of vectors for Sinkhorn divergence (primary only)
- `pe_diff_all`, `pe_diff_primary`: Vectors of PE differences
- `pe_ratio_all`, `pe_ratio_primary`: Vectors of PE ratios
- `prediction_success_all`, `prediction_success_primary`: Vectors of success flags

**MC Truth branches (only present when using `-mc` flag):**
- `vtx_dist_to_true`: Vector of distances from each vertex candidate to true neutrino vertex (cm)
- `true_vtx_x`, `true_vtx_y`, `true_vtx_z`: True neutrino vertex position (event-level)
- `has_mc_truth`: Boolean flag indicating if MC truth information was successfully extracted

## Using as Friend Tree

### Vectorized Format
```cpp
// In ROOT
TFile* f1 = TFile::Open("original_ntuple.root");
TTree* t1 = (TTree*)f1->Get("EventTree");
TFile* f2 = TFile::Open("flash_predictions_vectorized_friend.root");
TTree* t2 = (TTree*)f2->Get("FlashPredictionFriend");
t1->AddFriend(t2);

// Access first vertex candidate in each event
t1->Draw("pred_total_pe_all[0]:obs_total_pe", "n_vertices>0 && has_flash");

// Access all vertex candidates (requires loop or special syntax)
t1->Draw("pred_total_pe_all:obs_total_pe", "has_flash", "para");
```

### Original Format
```cpp
// In ROOT - requires more complex matching since entries don't align
// Generally recommend using vectorized format for friend trees
```

## Physics Parameters

The flash predictor uses the following parameters:
- ADC to electron conversion: 200 ADC/e⁻
- Energy per electron: 23.6 keV
- Photons per MeV: 24,000
- Recombination factor: 0.7
- Track step size: 0.3-0.5 cm
- Pixel window: 3×3 for charge collection

## Analysis Examples

### Vectorized Format
1. **PE ratio distribution for first vertex**:
   ```cpp
   tree->Draw("pe_ratio_all[0]", "n_vertices>0 && has_flash && obs_total_pe>50");
   ```

2. **Sinkhorn divergence vs PE difference for all vertices**:
   ```cpp
   tree->Draw("sinkhorn_div_all[0]:pe_diff_all", 
              "n_vertices>0 && has_flash", "para colz");
   ```

3. **Primary vs all particles comparison**:
   ```cpp
   tree->Draw("pred_total_pe_primary[0]:pred_total_pe_all[0]", 
              "n_vertices>0 && prediction_success_all[0] && prediction_success_primary[0]", "colz");
   ```

4. **Multi-vertex event analysis**:
   ```cpp
   tree->Draw("n_vertices", "has_vertices");
   tree->Draw("pred_total_pe_all", "n_vertices>1", "para"); // Events with multiple vertices
   ```

5. **MC Truth analysis (requires `-mc` flag during processing)**:
   ```cpp
   // Distance to true vertex vs Sinkhorn divergence
   tree->Draw("sinkhorn_div_all[0]:vtx_dist_to_true[0]", 
              "has_mc_truth && n_vertices>0 && prediction_success_all[0]", "colz");
   
   // Find best vertex candidate (closest to true)
   tree->Draw("vtx_dist_to_true", "has_mc_truth", "");
   
   // Correlation between distance to true vertex and PE ratio
   tree->Draw("pe_ratio_all[0]:vtx_dist_to_true[0]", 
              "has_mc_truth && n_vertices>0 && prediction_success_all[0]", "colz");
   ```

### Original Format
1. **PE ratio distribution**:
   ```cpp
   tree->Draw("pe_ratio_all", "has_flash && obs_total_pe>50");
   ```

2. **Sinkhorn divergence vs PE difference**:
   ```cpp
   tree->Draw("sinkhorn_div_all[1]:pe_diff_all", 
              "has_flash && prediction_success_all", "colz");
   ```

## TO DO:
1. ~~Make executable that calculates predictions and Sinkhorn divergence~~ ✓
2. ~~Write scripts for grid job submission~~ ✓
3. ~~Create combination script for friend tree creation~~ ✓
4. Perform analysis to optimize cuts on Sinkhorn divergence and PE difference