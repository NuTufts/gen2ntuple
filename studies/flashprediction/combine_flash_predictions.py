#!/usr/bin/env python3
"""
Script to combine vectorized flash prediction output files and create friend trees

This script:
1. Takes multiple vectorized flash prediction ROOT files
2. Combines them in the correct order to match the ntuple
3. Handles missing entries by filling with default values
4. Creates a friend tree that can be used with existing ntuples

The vectorized format stores all vertex candidates per event in vectors,
making it directly compatible with ntuple friend tree structure.
"""

import ROOT
import os
import sys
import argparse
import numpy as np
from collections import defaultdict

def get_entry_info(flash_file):
    """Get entry information from a vectorized flash prediction file"""
    
    tf = ROOT.TFile(flash_file, 'READ')
    tree = tf.Get('FlashPredictionTree')
    
    if not tree:
        print(f"Warning: No FlashPredictionTree found in {flash_file}")
        return {}
    
    entry_info = {}
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        entry_info[tree.entry] = {
            'run': tree.run,
            'subrun': tree.subrun,
            'event': tree.event,
            'n_vertices': tree.n_vertices
        }
    
    tf.Close()
    return entry_info

def create_combined_tree(output_file, input_files, ntuple_file=None):
    """
    Create combined flash prediction tree from vectorized files
    
    Args:
        output_file: Output ROOT file path
        input_files: List of input vectorized flash prediction files
        ntuple_file: Original ntuple file to match structure (optional)
    """
    
    print(f"Creating combined vectorized output: {output_file}")
    
    # First pass: determine which entries we have
    all_entries = set()
    entry_to_file = {}
    
    for input_file in input_files:
        print(f"Scanning {input_file}...")
        entry_info = get_entry_info(input_file)
        for entry in entry_info:
            all_entries.add(entry)
            entry_to_file[entry] = input_file
    
    print(f"Found flash predictions for {len(all_entries)} unique entries")
    
    # Determine total entries needed
    if ntuple_file:
        # Match ntuple structure
        ntuple_tf = ROOT.TFile(ntuple_file, 'READ')
        ntuple_tree = ntuple_tf.Get('EventTree')
        if ntuple_tree:
            total_entries = ntuple_tree.GetEntries()
            print(f"Matching ntuple with {total_entries} entries")
        else:
            total_entries = max(all_entries) + 1 if all_entries else 0
        ntuple_tf.Close()
    else:
        total_entries = max(all_entries) + 1 if all_entries else 0
    
    # Create output file and tree
    out_tf = ROOT.TFile(output_file, 'RECREATE')
    out_tree = ROOT.TTree('FlashPredictionFriend', 'Combined vectorized flash predictions for friend tree')
    
    # Define branches using same structure as vectorized version
    # Event-level info
    entry = np.zeros(1, dtype=np.int32)
    run = np.zeros(1, dtype=np.int32)
    subrun = np.zeros(1, dtype=np.int32)
    event = np.zeros(1, dtype=np.int32)
    n_vertices = np.zeros(1, dtype=np.int32)
    has_vertices = np.zeros(1, dtype=np.bool_)
    has_flash = np.zeros(1, dtype=np.bool_)
    
    # Observed flash info
    obs_total_pe = np.zeros(1, dtype=np.float32)
    obs_time = np.zeros(1, dtype=np.float32)
    
    # Set up branches for vectors (ROOT will handle the std::vector automatically)
    # We'll use the ROOT vector types directly
    obs_pe_per_pmt = ROOT.std.vector('float')()
    
    pred_total_pe_all_v = ROOT.std.vector('float')()
    pred_pe_per_pmt_all_v = ROOT.std.vector(ROOT.std.vector('float'))()
    n_tracks_all_v = ROOT.std.vector('int')()
    n_showers_all_v = ROOT.std.vector('int')()
    total_charge_all_v = ROOT.std.vector('float')()
    total_photons_all_v = ROOT.std.vector('float')()
    
    pred_total_pe_primary_v = ROOT.std.vector('float')()
    pred_pe_per_pmt_primary_v = ROOT.std.vector(ROOT.std.vector('float'))()
    n_tracks_primary_v = ROOT.std.vector('int')()
    n_showers_primary_v = ROOT.std.vector('int')()
    total_charge_primary_v = ROOT.std.vector('float')()
    total_photons_primary_v = ROOT.std.vector('float')()
    
    sinkhorn_div_all_v = ROOT.std.vector(ROOT.std.vector('float'))()
    sinkhorn_div_primary_v = ROOT.std.vector(ROOT.std.vector('float'))()
    pe_diff_all_v = ROOT.std.vector('float')()
    pe_diff_primary_v = ROOT.std.vector('float')()
    pe_ratio_all_v = ROOT.std.vector('float')()
    pe_ratio_primary_v = ROOT.std.vector('float')()
    
    prediction_success_all_v = ROOT.std.vector('bool')()
    prediction_success_primary_v = ROOT.std.vector('bool')()
    
    # Set up branches
    out_tree.Branch('entry', entry, 'entry/I')
    out_tree.Branch('run', run, 'run/I')
    out_tree.Branch('subrun', subrun, 'subrun/I')
    out_tree.Branch('event', event, 'event/I')
    out_tree.Branch('n_vertices', n_vertices, 'n_vertices/I')
    out_tree.Branch('has_vertices', has_vertices, 'has_vertices/O')
    out_tree.Branch('has_flash', has_flash, 'has_flash/O')
    
    out_tree.Branch('obs_total_pe', obs_total_pe, 'obs_total_pe/F')
    out_tree.Branch('obs_time', obs_time, 'obs_time/F')
    out_tree.Branch('obs_pe_per_pmt', obs_pe_per_pmt)
    
    out_tree.Branch('pred_total_pe_all', pred_total_pe_all_v)
    out_tree.Branch('pred_pe_per_pmt_all', pred_pe_per_pmt_all_v)
    out_tree.Branch('n_tracks_all', n_tracks_all_v)
    out_tree.Branch('n_showers_all', n_showers_all_v)
    out_tree.Branch('total_charge_all', total_charge_all_v)
    out_tree.Branch('total_photons_all', total_photons_all_v)
    
    out_tree.Branch('pred_total_pe_primary', pred_total_pe_primary_v)
    out_tree.Branch('pred_pe_per_pmt_primary', pred_pe_per_pmt_primary_v)
    out_tree.Branch('n_tracks_primary', n_tracks_primary_v)
    out_tree.Branch('n_showers_primary', n_showers_primary_v)
    out_tree.Branch('total_charge_primary', total_charge_primary_v)
    out_tree.Branch('total_photons_primary', total_photons_primary_v)
    
    out_tree.Branch('sinkhorn_div_all', sinkhorn_div_all_v)
    out_tree.Branch('sinkhorn_div_primary', sinkhorn_div_primary_v)
    out_tree.Branch('pe_diff_all', pe_diff_all_v)
    out_tree.Branch('pe_diff_primary', pe_diff_primary_v)
    out_tree.Branch('pe_ratio_all', pe_ratio_all_v)
    out_tree.Branch('pe_ratio_primary', pe_ratio_primary_v)
    
    out_tree.Branch('prediction_success_all', prediction_success_all_v)
    out_tree.Branch('prediction_success_primary', prediction_success_primary_v)
    
    # Process entries
    for ientry in range(total_entries):
        
        if ientry % 1000 == 0:
            print(f"Processing entry {ientry}/{total_entries}")
        
        # Clear vectors
        obs_pe_per_pmt.clear()
        pred_total_pe_all_v.clear()
        pred_pe_per_pmt_all_v.clear()
        n_tracks_all_v.clear()
        n_showers_all_v.clear()
        total_charge_all_v.clear()
        total_photons_all_v.clear()
        
        pred_total_pe_primary_v.clear()
        pred_pe_per_pmt_primary_v.clear()
        n_tracks_primary_v.clear()
        n_showers_primary_v.clear()
        total_charge_primary_v.clear()
        total_photons_primary_v.clear()
        
        sinkhorn_div_all_v.clear()
        sinkhorn_div_primary_v.clear()
        pe_diff_all_v.clear()
        pe_diff_primary_v.clear()
        pe_ratio_all_v.clear()
        pe_ratio_primary_v.clear()
        
        prediction_success_all_v.clear()
        prediction_success_primary_v.clear()
        
        # Set default values
        entry[0] = ientry
        run[0] = -1
        subrun[0] = -1
        event[0] = -1
        n_vertices[0] = 0
        has_vertices[0] = False
        has_flash[0] = False
        obs_total_pe[0] = -1
        obs_time[0] = -999
        
        # Fill obs_pe_per_pmt with default values
        for i in range(32):
            obs_pe_per_pmt.push_back(-1.0)
        
        # Check if we have data for this entry
        if ientry in entry_to_file:
            # Load data from the appropriate file
            input_file = entry_to_file[ientry]
            tf = ROOT.TFile(input_file, 'READ')
            tree = tf.Get('FlashPredictionTree')
            
            # Find the entry in the input tree
            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                
                if tree.entry != ientry:
                    continue
                
                # Copy event-level data
                run[0] = tree.run
                subrun[0] = tree.subrun
                event[0] = tree.event
                n_vertices[0] = tree.n_vertices
                has_vertices[0] = tree.has_vertices
                has_flash[0] = tree.has_flash
                obs_total_pe[0] = tree.obs_total_pe
                obs_time[0] = tree.obs_time
                
                # Copy observed PMT data
                obs_pe_per_pmt.clear()
                for j in range(tree.obs_pe_per_pmt.size()):
                    obs_pe_per_pmt.push_back(tree.obs_pe_per_pmt[j])
                
                # Copy vertex-level vectors directly
                for j in range(tree.pred_total_pe_all.size()):
                    pred_total_pe_all_v.push_back(tree.pred_total_pe_all[j])
                    
                for j in range(tree.pred_pe_per_pmt_all.size()):
                    pmt_vec = ROOT.std.vector('float')()
                    for k in range(tree.pred_pe_per_pmt_all[j].size()):
                        pmt_vec.push_back(tree.pred_pe_per_pmt_all[j][k])
                    pred_pe_per_pmt_all_v.push_back(pmt_vec)
                
                for j in range(tree.n_tracks_all.size()):
                    n_tracks_all_v.push_back(tree.n_tracks_all[j])
                    
                for j in range(tree.n_showers_all.size()):
                    n_showers_all_v.push_back(tree.n_showers_all[j])
                    
                for j in range(tree.total_charge_all.size()):
                    total_charge_all_v.push_back(tree.total_charge_all[j])
                    
                for j in range(tree.total_photons_all.size()):
                    total_photons_all_v.push_back(tree.total_photons_all[j])
                
                # Copy primary particle vectors
                for j in range(tree.pred_total_pe_primary.size()):
                    pred_total_pe_primary_v.push_back(tree.pred_total_pe_primary[j])
                    
                for j in range(tree.pred_pe_per_pmt_primary.size()):
                    pmt_vec = ROOT.std.vector('float')()
                    for k in range(tree.pred_pe_per_pmt_primary[j].size()):
                        pmt_vec.push_back(tree.pred_pe_per_pmt_primary[j][k])
                    pred_pe_per_pmt_primary_v.push_back(pmt_vec)
                
                for j in range(tree.n_tracks_primary.size()):
                    n_tracks_primary_v.push_back(tree.n_tracks_primary[j])
                    
                for j in range(tree.n_showers_primary.size()):
                    n_showers_primary_v.push_back(tree.n_showers_primary[j])
                    
                for j in range(tree.total_charge_primary.size()):
                    total_charge_primary_v.push_back(tree.total_charge_primary[j])
                    
                for j in range(tree.total_photons_primary.size()):
                    total_photons_primary_v.push_back(tree.total_photons_primary[j])
                
                # Copy metrics vectors
                for j in range(tree.sinkhorn_div_all.size()):
                    sinkhorn_vec = ROOT.std.vector('float')()
                    for k in range(tree.sinkhorn_div_all[j].size()):
                        sinkhorn_vec.push_back(tree.sinkhorn_div_all[j][k])
                    sinkhorn_div_all_v.push_back(sinkhorn_vec)
                    
                for j in range(tree.sinkhorn_div_primary.size()):
                    sinkhorn_vec = ROOT.std.vector('float')()
                    for k in range(tree.sinkhorn_div_primary[j].size()):
                        sinkhorn_vec.push_back(tree.sinkhorn_div_primary[j][k])
                    sinkhorn_div_primary_v.push_back(sinkhorn_vec)
                
                for j in range(tree.pe_diff_all.size()):
                    pe_diff_all_v.push_back(tree.pe_diff_all[j])
                    
                for j in range(tree.pe_diff_primary.size()):
                    pe_diff_primary_v.push_back(tree.pe_diff_primary[j])
                    
                for j in range(tree.pe_ratio_all.size()):
                    pe_ratio_all_v.push_back(tree.pe_ratio_all[j])
                    
                for j in range(tree.pe_ratio_primary.size()):
                    pe_ratio_primary_v.push_back(tree.pe_ratio_primary[j])
                
                for j in range(tree.prediction_success_all.size()):
                    prediction_success_all_v.push_back(tree.prediction_success_all[j])
                    
                for j in range(tree.prediction_success_primary.size()):
                    prediction_success_primary_v.push_back(tree.prediction_success_primary[j])
                
                break
            
            tf.Close()
        
        # Fill tree
        out_tree.Fill()
    
    # Write and close
    out_tf.cd()
    out_tree.Write()
    
    print(f"\nCombined tree summary:")
    print(f"Total entries: {out_tree.GetEntries()}")
    print(f"Output file: {output_file}")
    
    out_tf.Close()

def main():
    parser = argparse.ArgumentParser(description="Combine vectorized flash prediction outputs")
    
    parser.add_argument("-i", "--input-files", nargs='+', required=True,
                        help="Input vectorized flash prediction ROOT files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output combined ROOT file")
    parser.add_argument("-n", "--ntuple", default=None,
                        help="Original ntuple file to match structure")
    
    args = parser.parse_args()
    
    # Check input files exist
    for f in args.input_files:
        if not os.path.exists(f):
            print(f"Error: Input file not found: {f}")
            sys.exit(1)
    
    # Create combined tree
    create_combined_tree(args.output, args.input_files, args.ntuple)
    
    print("\nTo use as a friend tree:")
    print("  TFile* f1 = TFile::Open('original_ntuple.root');")
    print("  TTree* t1 = (TTree*)f1->Get('EventTree');")
    print(f"  TFile* f2 = TFile::Open('{args.output}');")
    print("  TTree* t2 = (TTree*)f2->Get('FlashPredictionFriend');")
    print("  t1->AddFriend(t2);")
    print("\nExample analysis:")
    print("  t1->Draw('pred_total_pe_all[0]:obs_total_pe', 'n_vertices>0 && has_flash');")

if __name__ == "__main__":
    main()