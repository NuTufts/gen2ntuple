#!/usr/bin/env python3
"""
Quick summary script for flash prediction ROOT files
Shows a compact table of the data
"""

import ROOT
import sys
import argparse

def print_compact_summary(filename, max_entries=10):
    """Print a compact summary table of the flash prediction data"""
    
    # Open file and get tree
    f = ROOT.TFile.Open(filename, "READ")
    if not f or f.IsZombie():
        print(f"Error: Cannot open {filename}")
        return
        
    tree = f.Get("FlashPredictionTree")
    if not tree:
        print("Error: No FlashPredictionTree found")
        return
    
    print(f"\nFlash Prediction Summary for: {filename}")
    print(f"Total entries: {tree.GetEntries()}")
    print("\n" + "="*120)
    print(f"{'Entry':>5} {'Run':>6} {'SubRun':>6} {'Event':>6} {'#Vtx':>4} "
          f"{'ObsPE':>8} {'Vtx0_PE':>8} {'Ratio0':>7} {'Sinkhorn0':>10} "
          f"{'Vtx1_PE':>8} {'Ratio1':>7} {'Status':>10}")
    print("="*120)
    
    # Loop through entries
    n_entries = min(max_entries, tree.GetEntries())
    
    for i in range(n_entries):
        tree.GetEntry(i)
        
        # Basic info
        entry_str = f"{tree.entry:5d} {tree.run:6d} {tree.subrun:6d} {tree.event:6d} {tree.n_vertices:4d}"
        
        # Observed flash
        if tree.has_flash:
            obs_pe_str = f"{tree.obs_total_pe:8.1f}"
        else:
            obs_pe_str = f"{'--':>8}"
        
        # First vertex
        if tree.n_vertices > 0 and len(tree.pred_total_pe_all) > 0:
            vtx0_pe = tree.pred_total_pe_all[0]
            vtx0_pe_str = f"{vtx0_pe:8.1f}"
            
            if tree.has_flash and tree.obs_total_pe > 0:
                ratio0 = vtx0_pe / tree.obs_total_pe
                ratio0_str = f"{ratio0:7.3f}"
            else:
                ratio0_str = f"{'--':>7}"
            
            # Sinkhorn divergence (λ=1.0)
            if len(tree.sinkhorn_div_all) > 0 and len(tree.sinkhorn_div_all[0]) > 1:
                sinkhorn0 = tree.sinkhorn_div_all[0][1]  # λ=1.0
                if sinkhorn0 > 0:
                    sinkhorn0_str = f"{sinkhorn0:10.2f}"
                else:
                    sinkhorn0_str = f"{'--':>10}"
            else:
                sinkhorn0_str = f"{'--':>10}"
        else:
            vtx0_pe_str = f"{'--':>8}"
            ratio0_str = f"{'--':>7}"
            sinkhorn0_str = f"{'--':>10}"
        
        # Second vertex (if exists)
        if tree.n_vertices > 1 and len(tree.pred_total_pe_all) > 1:
            vtx1_pe = tree.pred_total_pe_all[1]
            vtx1_pe_str = f"{vtx1_pe:8.1f}"
            
            if tree.has_flash and tree.obs_total_pe > 0:
                ratio1 = vtx1_pe / tree.obs_total_pe
                ratio1_str = f"{ratio1:7.3f}"
            else:
                ratio1_str = f"{'--':>7}"
        else:
            vtx1_pe_str = f"{'--':>8}"
            ratio1_str = f"{'--':>7}"
        
        # Status
        if tree.has_vertices and tree.has_flash:
            status = "Both"
        elif tree.has_vertices:
            status = "VtxOnly"
        elif tree.has_flash:
            status = "FlashOnly"
        else:
            status = "Neither"
        
        # Print row
        print(f"{entry_str} {obs_pe_str} {vtx0_pe_str} {ratio0_str} {sinkhorn0_str} "
              f"{vtx1_pe_str} {ratio1_str} {status:>10}")
    
    if tree.GetEntries() > max_entries:
        print(f"... ({tree.GetEntries() - max_entries} more entries)")
    
    print("="*120)
    
    # Quick statistics
    total_vertices = 0
    successful_predictions = 0
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        total_vertices += tree.n_vertices
        successful_predictions += sum(tree.prediction_success_all)
    
    print(f"\nQuick Stats:")
    print(f"  Total vertices: {total_vertices}")
    print(f"  Successful predictions: {successful_predictions} ({100*successful_predictions/max(1,total_vertices):.1f}%)")
    print(f"  Average vertices/event: {total_vertices/tree.GetEntries():.2f}")
    
    f.Close()

def main():
    parser = argparse.ArgumentParser(description="Quick summary of flash prediction ROOT file")
    parser.add_argument("input_file", help="Input ROOT file")
    parser.add_argument("-n", "--num-entries", type=int, default=20,
                        help="Maximum entries to show (default: 20)")
    
    args = parser.parse_args()
    
    print_compact_summary(args.input_file, args.num_entries)

if __name__ == "__main__":
    main()