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
  [-v]
```

**Parameters:**
- `-d, --dlmerged`: Input dlmerged ROOT file (contains ADC images and opflash)
- `-r, --reco`: Input reco file with KPSRecoManagerTree
- `-o, --output`: Output ROOT file for flash predictions
- `-n, --num-entries`: Number of entries to process (default: all)
- `-s, --start-entry`: Starting entry (default: 0)
- `-t, --threshold`: ADC threshold (default: 10.0)
- `-v, --verbose`: Enable verbose output

### Submitting batch jobs

For processing many files on a cluster:

```bash
python3 submit_flash_prediction_jobs.py \
  -d "/path/to/dlmerged*.root" \
  -r "/path/to/reco*.root" \
  -o /path/to/output/dir \
  [-j jobs_per_file] \
  [-e entries_per_job] \
  [--dry-run]
```

### Combining outputs

To create a single friend tree from multiple output files:

```bash
python3 combine_flash_predictions.py \
  -i output1.root output2.root ... \
  -o combined_friend_tree.root \
  [-n original_ntuple.root] \
  [-m max_vertices_per_entry]
```

## Output Structure

The output ROOT tree contains the following branches per vertex:

**Predictions:**
- `pred_total_pe_all`: Total predicted PE using all particles
- `pred_total_pe_primary`: Total predicted PE using primary particles only
- `pred_pe_per_pmt_all[32]`: Per-PMT predictions (all particles)
- `pred_pe_per_pmt_primary[32]`: Per-PMT predictions (primary only)

**Particle counts:**
- `n_tracks_all`, `n_showers_all`: Track/shower count (all particles)
- `n_tracks_primary`, `n_showers_primary`: Track/shower count (primary only)

**Metrics:**
- `sinkhorn_div_all[3]`: Sinkhorn divergence for λ = 0.1, 1.0, 10.0 (all)
- `sinkhorn_div_primary[3]`: Sinkhorn divergence (primary only)
- `pe_diff_all/primary`: Predicted - Observed total PE
- `pe_ratio_all/primary`: Predicted / Observed total PE

**Observed flash:**
- `obs_total_pe`: Observed total PE
- `obs_pe_per_pmt[32]`: Observed PE per PMT
- `obs_time`: Flash time

**Status flags:**
- `has_vertex`: Whether vertex candidates exist
- `has_flash`: Whether observed flash exists
- `prediction_success_all/primary`: Whether prediction succeeded

## Using as Friend Tree

```cpp
// In ROOT
TFile* f1 = TFile::Open("original_ntuple.root");
TTree* t1 = (TTree*)f1->Get("EventTree");
TFile* f2 = TFile::Open("flash_predictions_friend.root");
TTree* t2 = (TTree*)f2->Get("FlashPredictionFriend");
t1->AddFriend(t2);

// Now you can access flash variables in your analysis
t1->Draw("pred_total_pe_all:obs_total_pe", "has_flash && prediction_success_all");
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

1. **PE ratio distribution**:
   ```cpp
   tree->Draw("pe_ratio_all", "has_flash && obs_total_pe>50");
   ```

2. **Sinkhorn divergence vs PE difference**:
   ```cpp
   tree->Draw("sinkhorn_div_all[1]:pe_diff_all", 
              "has_flash && prediction_success_all", "colz");
   ```

3. **Primary vs all particles comparison**:
   ```cpp
   tree->Draw("pred_total_pe_primary:pred_total_pe_all", 
              "prediction_success_all && prediction_success_primary", "colz");
   ```

## TO DO:
1. ~~Make executable that calculates predictions and Sinkhorn divergence~~ ✓
2. ~~Write scripts for grid job submission~~ ✓
3. ~~Create combination script for friend tree creation~~ ✓
4. Perform analysis to optimize cuts on Sinkhorn divergence and PE difference