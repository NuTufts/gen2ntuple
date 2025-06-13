#!/usr/bin/env python3
"""
Script to combine flash prediction output files and create friend trees

This script:
1. Takes multiple flash prediction ROOT files
2. Combines them in the correct order to match the ntuple
3. Handles missing entries by filling with default values
4. Creates a friend tree that can be used with existing ntuples
"""

import ROOT
import os
import sys
import argparse
import numpy as np
from collections import defaultdict

def get_entry_mapping(flash_file):
    """Get mapping of (entry, vertex_idx) from a flash prediction file"""
    
    tf = ROOT.TFile(flash_file, 'READ')
    tree = tf.Get('FlashPredictionTree')
    
    if not tree:
        print(f"Warning: No FlashPredictionTree found in {flash_file}")
        return {}
    
    mapping = defaultdict(list)
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        mapping[tree.entry].append(tree.vertex_idx)
    
    tf.Close()
    return mapping

def create_combined_tree(output_file, input_files, ntuple_file=None, max_vertices=10):
    """
    Create combined flash prediction tree
    
    Args:
        output_file: Output ROOT file path
        input_files: List of input flash prediction files
        ntuple_file: Original ntuple file to match structure (optional)
        max_vertices: Maximum vertices per entry to store
    """
    
    print(f"Creating combined output: {output_file}")
    
    # First pass: determine which entries we have
    all_entries = set()
    entry_to_file = {}
    
    for input_file in input_files:
        print(f"Scanning {input_file}...")
        mapping = get_entry_mapping(input_file)
        for entry in mapping:
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
    out_tree = ROOT.TTree('FlashPredictionFriend', 'Combined flash predictions for friend tree')
    
    # Define arrays for storing multiple vertices per entry
    # Using std::vector would be cleaner but numpy arrays work for now
    n_vertices = np.zeros(1, dtype=np.int32)
    
    # Arrays for each vertex (up to max_vertices)
    vertex_idx = np.zeros(max_vertices, dtype=np.int32)
    pred_total_pe_all = np.zeros(max_vertices, dtype=np.float32)
    pred_total_pe_primary = np.zeros(max_vertices, dtype=np.float32)
    n_tracks_all = np.zeros(max_vertices, dtype=np.int32)
    n_showers_all = np.zeros(max_vertices, dtype=np.int32)
    n_tracks_primary = np.zeros(max_vertices, dtype=np.int32)
    n_showers_primary = np.zeros(max_vertices, dtype=np.int32)
    
    # Sinkhorn divergences (3 regularization values per vertex)
    sinkhorn_div_all = np.zeros((max_vertices, 3), dtype=np.float32)
    sinkhorn_div_primary = np.zeros((max_vertices, 3), dtype=np.float32)
    
    # PE differences and ratios
    pe_diff_all = np.zeros(max_vertices, dtype=np.float32)
    pe_diff_primary = np.zeros(max_vertices, dtype=np.float32)
    pe_ratio_all = np.zeros(max_vertices, dtype=np.float32)
    pe_ratio_primary = np.zeros(max_vertices, dtype=np.float32)
    
    # Status flags
    has_vertex = np.zeros(max_vertices, dtype=np.bool_)
    has_flash = np.zeros(max_vertices, dtype=np.bool_)
    prediction_success_all = np.zeros(max_vertices, dtype=np.bool_)
    prediction_success_primary = np.zeros(max_vertices, dtype=np.bool_)
    
    # Observed flash info (same for all vertices in an entry)
    obs_total_pe = np.zeros(1, dtype=np.float32)
    obs_time = np.zeros(1, dtype=np.float32)
    
    # Set up branches
    out_tree.Branch('n_vertices', n_vertices, 'n_vertices/I')
    out_tree.Branch('vertex_idx', vertex_idx, 'vertex_idx[n_vertices]/I')
    
    out_tree.Branch('pred_total_pe_all', pred_total_pe_all, 'pred_total_pe_all[n_vertices]/F')
    out_tree.Branch('pred_total_pe_primary', pred_total_pe_primary, 'pred_total_pe_primary[n_vertices]/F')
    
    out_tree.Branch('n_tracks_all', n_tracks_all, 'n_tracks_all[n_vertices]/I')
    out_tree.Branch('n_showers_all', n_showers_all, 'n_showers_all[n_vertices]/I')
    out_tree.Branch('n_tracks_primary', n_tracks_primary, 'n_tracks_primary[n_vertices]/I')
    out_tree.Branch('n_showers_primary', n_showers_primary, 'n_showers_primary[n_vertices]/I')
    
    # Sinkhorn divergences (flattened)
    out_tree.Branch('sinkhorn_div_all_0p1', sinkhorn_div_all[:, 0], 'sinkhorn_div_all_0p1[n_vertices]/F')
    out_tree.Branch('sinkhorn_div_all_1p0', sinkhorn_div_all[:, 1], 'sinkhorn_div_all_1p0[n_vertices]/F')
    out_tree.Branch('sinkhorn_div_all_10p0', sinkhorn_div_all[:, 2], 'sinkhorn_div_all_10p0[n_vertices]/F')
    
    out_tree.Branch('sinkhorn_div_primary_0p1', sinkhorn_div_primary[:, 0], 'sinkhorn_div_primary_0p1[n_vertices]/F')
    out_tree.Branch('sinkhorn_div_primary_1p0', sinkhorn_div_primary[:, 1], 'sinkhorn_div_primary_1p0[n_vertices]/F')
    out_tree.Branch('sinkhorn_div_primary_10p0', sinkhorn_div_primary[:, 2], 'sinkhorn_div_primary_10p0[n_vertices]/F')
    
    out_tree.Branch('pe_diff_all', pe_diff_all, 'pe_diff_all[n_vertices]/F')
    out_tree.Branch('pe_diff_primary', pe_diff_primary, 'pe_diff_primary[n_vertices]/F')
    out_tree.Branch('pe_ratio_all', pe_ratio_all, 'pe_ratio_all[n_vertices]/F')
    out_tree.Branch('pe_ratio_primary', pe_ratio_primary, 'pe_ratio_primary[n_vertices]/F')
    
    out_tree.Branch('has_vertex', has_vertex, 'has_vertex[n_vertices]/O')
    out_tree.Branch('has_flash', has_flash, 'has_flash[n_vertices]/O')
    out_tree.Branch('prediction_success_all', prediction_success_all, 'prediction_success_all[n_vertices]/O')
    out_tree.Branch('prediction_success_primary', prediction_success_primary, 'prediction_success_primary[n_vertices]/O')
    
    out_tree.Branch('obs_total_pe', obs_total_pe, 'obs_total_pe/F')
    out_tree.Branch('obs_time', obs_time, 'obs_time/F')
    
    # Process entries
    for entry in range(total_entries):
        
        if entry % 1000 == 0:
            print(f"Processing entry {entry}/{total_entries}")
        
        # Reset arrays
        n_vertices[0] = 0
        obs_total_pe[0] = -1
        obs_time[0] = -999
        
        # Clear vertex arrays
        for i in range(max_vertices):
            vertex_idx[i] = -1
            pred_total_pe_all[i] = -1
            pred_total_pe_primary[i] = -1
            n_tracks_all[i] = 0
            n_showers_all[i] = 0
            n_tracks_primary[i] = 0
            n_showers_primary[i] = 0
            pe_diff_all[i] = -999
            pe_diff_primary[i] = -999
            pe_ratio_all[i] = -999
            pe_ratio_primary[i] = -999
            has_vertex[i] = False
            has_flash[i] = False
            prediction_success_all[i] = False
            prediction_success_primary[i] = False
            
            for j in range(3):
                sinkhorn_div_all[i, j] = -999
                sinkhorn_div_primary[i, j] = -999
        
        # Check if we have data for this entry
        if entry in entry_to_file:
            # Load data from the appropriate file
            input_file = entry_to_file[entry]
            tf = ROOT.TFile(input_file, 'READ')
            tree = tf.Get('FlashPredictionTree')
            
            vtx_count = 0
            
            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                
                if tree.entry != entry:
                    continue
                
                if vtx_count >= max_vertices:
                    print(f"Warning: Entry {entry} has more than {max_vertices} vertices, truncating")
                    break
                
                # Copy data
                vertex_idx[vtx_count] = tree.vertex_idx
                pred_total_pe_all[vtx_count] = tree.pred_total_pe_all
                pred_total_pe_primary[vtx_count] = tree.pred_total_pe_primary
                n_tracks_all[vtx_count] = tree.n_tracks_all
                n_showers_all[vtx_count] = tree.n_showers_all
                n_tracks_primary[vtx_count] = tree.n_tracks_primary
                n_showers_primary[vtx_count] = tree.n_showers_primary
                
                # Copy sinkhorn divergences
                for j in range(3):
                    sinkhorn_div_all[vtx_count, j] = tree.sinkhorn_div_all[j]
                    sinkhorn_div_primary[vtx_count, j] = tree.sinkhorn_div_primary[j]
                
                pe_diff_all[vtx_count] = tree.pe_diff_all
                pe_diff_primary[vtx_count] = tree.pe_diff_primary
                pe_ratio_all[vtx_count] = tree.pe_ratio_all
                pe_ratio_primary[vtx_count] = tree.pe_ratio_primary
                
                has_vertex[vtx_count] = tree.has_vertex
                has_flash[vtx_count] = tree.has_flash
                prediction_success_all[vtx_count] = tree.prediction_success_all
                prediction_success_primary[vtx_count] = tree.prediction_success_primary
                
                # Observed flash info (same for all vertices)
                if vtx_count == 0:
                    obs_total_pe[0] = tree.obs_total_pe
                    obs_time[0] = tree.obs_time
                
                vtx_count += 1
            
            n_vertices[0] = vtx_count
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
    parser = argparse.ArgumentParser(description="Combine flash prediction outputs")
    
    parser.add_argument("-i", "--input-files", nargs='+', required=True,
                        help="Input flash prediction ROOT files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output combined ROOT file")
    parser.add_argument("-n", "--ntuple", default=None,
                        help="Original ntuple file to match structure")
    parser.add_argument("-m", "--max-vertices", type=int, default=10,
                        help="Maximum vertices per entry (default: 10)")
    
    args = parser.parse_args()
    
    # Check input files exist
    for f in args.input_files:
        if not os.path.exists(f):
            print(f"Error: Input file not found: {f}")
            sys.exit(1)
    
    # Create combined tree
    create_combined_tree(args.output, args.input_files, args.ntuple, args.max_vertices)
    
    print("\nTo use as a friend tree:")
    print("  TFile* f1 = TFile::Open('original_ntuple.root');")
    print("  TTree* t1 = (TTree*)f1->Get('EventTree');")
    print(f"  TFile* f2 = TFile::Open('{args.output}');")
    print("  TTree* t2 = (TTree*)f2->Get('FlashPredictionFriend');")
    print("  t1->AddFriend(t2);")

if __name__ == "__main__":
    main()