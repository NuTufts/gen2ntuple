#!/usr/bin/env python3
"""
Script to dump and analyze the contents of ROOT files created by calculate_flash_predictions_vectorized

This script will:
1. Open the ROOT file and access the FlashPredictionTree
2. Print the tree structure and branch information
3. Loop through entries and display the content
4. Provide summary statistics
"""

import ROOT
import sys
import argparse
import numpy as np
from collections import defaultdict

def print_branch_info(tree):
    """Print information about all branches in the tree"""
    print("\n" + "="*80)
    print("TREE STRUCTURE")
    print("="*80)
    print(f"Tree name: {tree.GetName()}")
    print(f"Tree title: {tree.GetTitle()}")
    print(f"Number of entries: {tree.GetEntries()}")
    
    print("\nBRANCHES:")
    print("-"*80)
    
    branches = tree.GetListOfBranches()
    for branch in branches:
        branch_name = branch.GetName()
        branch_class = branch.GetClassName()
        
        # Get leaf information
        leaves = branch.GetListOfLeaves()
        if leaves.GetEntries() > 0:
            leaf = leaves[0]
            leaf_type = leaf.GetTypeName()
            
            # Check if it's an array
            if leaf.GetLen() > 1 and not branch_class:
                print(f"{branch_name:<30} {leaf_type}[{leaf.GetLen()}]")
            elif branch_class:
                print(f"{branch_name:<30} {branch_class}")
            else:
                print(f"{branch_name:<30} {leaf_type}")
        else:
            print(f"{branch_name:<30} {branch_class}")

def dump_entry_details(tree, entry_num, max_vertices=3):
    """Dump detailed information for a specific entry"""
    tree.GetEntry(entry_num)
    
    print(f"\n{'='*80}")
    print(f"ENTRY {entry_num} DETAILS")
    print(f"{'='*80}")
    
    # Event identification
    print(f"\nEvent Info:")
    print(f"  Entry: {tree.entry}")
    print(f"  Run: {tree.run}, Subrun: {tree.subrun}, Event: {tree.event}")
    print(f"  Number of vertices: {tree.n_vertices}")
    print(f"  Has vertices: {tree.has_vertices}")
    print(f"  Has flash: {tree.has_flash}")
    
    # MC truth info (if available)
    if hasattr(tree, 'has_mc_truth'):
        print(f"  Has MC truth: {tree.has_mc_truth}")
        if tree.has_mc_truth:
            print(f"  True vertex: ({tree.true_vtx_x:.2f}, {tree.true_vtx_y:.2f}, {tree.true_vtx_z:.2f}) cm")
    
    # Observed flash info
    print(f"\nObserved Flash:")
    print(f"  Total PE: {tree.obs_total_pe:.2f}")
    print(f"  Time: {tree.obs_time:.2f} μs")
    
    # PMT distribution for observed flash
    if hasattr(tree, 'obs_pe_per_pmt') and tree.obs_pe_per_pmt.size() > 0:
        print(f"  PMT PE distribution (first 10 PMTs):")
        for pmt in range(min(10, tree.obs_pe_per_pmt.size())):
            print(f"    PMT {pmt:2d}: {tree.obs_pe_per_pmt[pmt]:6.2f} PE")
        if tree.obs_pe_per_pmt.size() > 10:
            print(f"    ... ({tree.obs_pe_per_pmt.size() - 10} more PMTs)")
    
    # Vertex-specific information
    if tree.n_vertices > 0:
        print(f"\nVertex Predictions:")
        
        # Limit vertices shown
        vertices_to_show = min(tree.n_vertices, max_vertices)
        
        for vtx in range(vertices_to_show):
            print(f"\n  Vertex {vtx}:")
            
            # All particles prediction
            if vtx < len(tree.pred_total_pe_all):
                print(f"    All particles:")
                print(f"      Total PE: {tree.pred_total_pe_all[vtx]:.2f}")
                print(f"      Tracks: {tree.n_tracks_all[vtx]}, Showers: {tree.n_showers_all[vtx]}")
                print(f"      Total charge: {tree.total_charge_all[vtx]:.1f} ADC")
                print(f"      Total photons: {tree.total_photons_all[vtx]:.0f}")
                print(f"      Success: {tree.prediction_success_all[vtx]}")
                
                # Metrics
                if tree.pe_diff_all[vtx] > -900:  # Check for valid value
                    print(f"      PE difference: {tree.pe_diff_all[vtx]:.2f}")
                    print(f"      PE ratio: {tree.pe_ratio_all[vtx]:.4f}")
                
                # Sinkhorn divergences
                if vtx < len(tree.sinkhorn_div_all):
                    print(f"      Sinkhorn divergences:")
                    for i, reg in enumerate([0.1, 1.0, 10.0]):
                        if i < len(tree.sinkhorn_div_all[vtx]):
                            div_val = tree.sinkhorn_div_all[vtx][i]
                            if div_val > -1:  # Valid value
                                print(f"        λ={reg:4.1f}: {div_val:8.4f}")
            
            # Primary particles prediction
            if vtx < len(tree.pred_total_pe_primary):
                print(f"    Primary particles only:")
                print(f"      Total PE: {tree.pred_total_pe_primary[vtx]:.2f}")
                print(f"      Tracks: {tree.n_tracks_primary[vtx]}, Showers: {tree.n_showers_primary[vtx]}")
                print(f"      Success: {tree.prediction_success_primary[vtx]}")
                
                if tree.pe_diff_primary[vtx] > -900:  # Check for valid value
                    print(f"      PE difference: {tree.pe_diff_primary[vtx]:.2f}")
                    print(f"      PE ratio: {tree.pe_ratio_primary[vtx]:.4f}")
            
            # MC truth distance (if available)
            if hasattr(tree, 'vtx_dist_to_true') and hasattr(tree, 'has_mc_truth'):
                if tree.has_mc_truth and vtx < len(tree.vtx_dist_to_true):
                    print(f"    Distance to true vertex: {tree.vtx_dist_to_true[vtx]:.2f} cm")
        
        if tree.n_vertices > max_vertices:
            print(f"\n  ... and {tree.n_vertices - max_vertices} more vertices")

def analyze_statistics(tree):
    """Compute and print statistics across all entries"""
    print(f"\n{'='*80}")
    print("STATISTICS ACROSS ALL ENTRIES")
    print(f"{'='*80}")
    
    # Counters
    total_vertices = 0
    events_with_vertices = 0
    events_with_flash = 0
    successful_predictions_all = 0
    successful_predictions_primary = 0
    
    # For statistics
    all_pred_pe_all = []
    all_pred_pe_primary = []
    all_pe_ratios_all = []
    all_sinkhorn_divs = defaultdict(list)
    vertices_per_event = []
    
    # Loop through all entries
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        
        if tree.has_vertices:
            events_with_vertices += 1
            vertices_per_event.append(tree.n_vertices)
        
        if tree.has_flash:
            events_with_flash += 1
        
        total_vertices += tree.n_vertices
        
        # Collect per-vertex statistics
        for vtx in range(tree.n_vertices):
            if vtx < len(tree.prediction_success_all) and tree.prediction_success_all[vtx]:
                successful_predictions_all += 1
                all_pred_pe_all.append(tree.pred_total_pe_all[vtx])
                
                if tree.obs_total_pe > 0:
                    all_pe_ratios_all.append(tree.pe_ratio_all[vtx])
                
                # Sinkhorn divergences
                if vtx < len(tree.sinkhorn_div_all):
                    for i in range(min(3, len(tree.sinkhorn_div_all[vtx]))):
                        if tree.sinkhorn_div_all[vtx][i] > -1:
                            all_sinkhorn_divs[i].append(tree.sinkhorn_div_all[vtx][i])
            
            if vtx < len(tree.prediction_success_primary) and tree.prediction_success_primary[vtx]:
                successful_predictions_primary += 1
                all_pred_pe_primary.append(tree.pred_total_pe_primary[vtx])
    
    # Print summary
    print(f"\nEvent Summary:")
    print(f"  Total entries: {tree.GetEntries()}")
    print(f"  Events with vertices: {events_with_vertices} ({100*events_with_vertices/tree.GetEntries():.1f}%)")
    print(f"  Events with flash: {events_with_flash} ({100*events_with_flash/tree.GetEntries():.1f}%)")
    
    print(f"\nVertex Summary:")
    print(f"  Total vertices: {total_vertices}")
    if vertices_per_event:
        print(f"  Vertices per event: min={min(vertices_per_event)}, max={max(vertices_per_event)}, "
              f"mean={np.mean(vertices_per_event):.2f}")
    
    print(f"\nPrediction Success:")
    print(f"  All particles: {successful_predictions_all}/{total_vertices} "
          f"({100*successful_predictions_all/max(1,total_vertices):.1f}%)")
    print(f"  Primary only: {successful_predictions_primary}/{total_vertices} "
          f"({100*successful_predictions_primary/max(1,total_vertices):.1f}%)")
    
    if all_pred_pe_all:
        print(f"\nPredicted PE Statistics (all particles):")
        print(f"  Min: {min(all_pred_pe_all):.2f}")
        print(f"  Max: {max(all_pred_pe_all):.2f}")
        print(f"  Mean: {np.mean(all_pred_pe_all):.2f}")
        print(f"  Median: {np.median(all_pred_pe_all):.2f}")
        print(f"  Std: {np.std(all_pred_pe_all):.2f}")
    
    if all_pe_ratios_all:
        print(f"\nPE Ratio Statistics (predicted/observed):")
        print(f"  Min: {min(all_pe_ratios_all):.4f}")
        print(f"  Max: {max(all_pe_ratios_all):.4f}")
        print(f"  Mean: {np.mean(all_pe_ratios_all):.4f}")
        print(f"  Median: {np.median(all_pe_ratios_all):.4f}")
    
    if all_sinkhorn_divs:
        print(f"\nSinkhorn Divergence Statistics:")
        reg_params = [0.1, 1.0, 10.0]
        for i, divs in all_sinkhorn_divs.items():
            if divs and i < len(reg_params):
                print(f"  λ={reg_params[i]}:")
                print(f"    Min: {min(divs):.4f}")
                print(f"    Max: {max(divs):.4f}")
                print(f"    Mean: {np.mean(divs):.4f}")
                print(f"    Median: {np.median(divs):.4f}")

def main():
    parser = argparse.ArgumentParser(description="Dump flash prediction ROOT file contents")
    parser.add_argument("input_file", help="Input ROOT file")
    parser.add_argument("-n", "--num-entries", type=int, default=5,
                        help="Number of entries to dump in detail (default: 5)")
    parser.add_argument("-s", "--start-entry", type=int, default=0,
                        help="Starting entry for detailed dump (default: 0)")
    parser.add_argument("-m", "--max-vertices", type=int, default=3,
                        help="Maximum vertices to show per entry (default: 3)")
    parser.add_argument("--stats-only", action="store_true",
                        help="Only show statistics, skip detailed dump")
    parser.add_argument("--no-stats", action="store_true",
                        help="Skip statistics calculation")
    
    args = parser.parse_args()
    
    # Open ROOT file
    print(f"Opening ROOT file: {args.input_file}")
    root_file = ROOT.TFile.Open(args.input_file, "READ")
    
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open file {args.input_file}")
        sys.exit(1)
    
    # Get tree
    tree = root_file.Get("FlashPredictionTree")
    if not tree:
        print("Error: Cannot find FlashPredictionTree in file")
        print("Available objects in file:")
        root_file.ls()
        sys.exit(1)
    
    # Print branch information
    print_branch_info(tree)
    
    # Dump detailed entry information
    if not args.stats_only:
        print(f"\n{'='*80}")
        print(f"DETAILED ENTRY DUMP (showing {args.num_entries} entries starting from {args.start_entry})")
        print(f"{'='*80}")
        
        end_entry = min(args.start_entry + args.num_entries, tree.GetEntries())
        for entry in range(args.start_entry, end_entry):
            dump_entry_details(tree, entry, args.max_vertices)
    
    # Calculate and print statistics
    if not args.no_stats:
        analyze_statistics(tree)
    
    root_file.Close()
    print("\nDone!")

if __name__ == "__main__":
    main()