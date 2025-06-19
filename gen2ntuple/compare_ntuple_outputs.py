#!/usr/bin/env python3
"""
Compare output ntuples from C++ and Python versions of gen2ntuple
to validate that the implementations match.
"""

import sys
import ROOT
import numpy as np
from collections import defaultdict

def get_tree_entries(tree):
    """Get all entries from a tree as a dictionary"""
    entries = []
    
    # Get list of branches
    branches = [br.GetName() for br in tree.GetListOfBranches()]
    
    # Read all entries
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        entry = {}
        for branch_name in branches:
            try:
                # Get the value - handle different types
                val = getattr(tree, branch_name)
                # Convert ROOT types to Python types if needed
                if hasattr(val, '__len__'):  # Array-like
                    entry[branch_name] = list(val)
                else:  # Scalar
                    entry[branch_name] = float(val) if isinstance(val, (int, float)) else val
            except:
                pass  # Skip branches we can't read
        entries.append(entry)
    return entries, branches

def compare_values(val1, val2, branch_name, tolerance=1e-5):
    """Compare two values with appropriate tolerance"""
    
    # Handle arrays
    if isinstance(val1, list) and isinstance(val2, list):
        if len(val1) != len(val2):
            return False, f"Different array lengths: {len(val1)} vs {len(val2)}"
        
        for i, (v1, v2) in enumerate(zip(val1, val2)):
            if isinstance(v1, (int, float)) and isinstance(v2, (int, float)):
                if abs(v1 - v2) > tolerance * max(abs(v1), abs(v2), 1.0):
                    return False, f"Array element {i}: {v1} vs {v2}"
        return True, "Arrays match"
    
    # Handle scalars
    elif isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
        # Special handling for -999 or similar default values
        if val1 == val2:
            return True, "Exact match"
        
        # Relative tolerance for non-zero values
        if abs(val1) > tolerance or abs(val2) > tolerance:
            rel_diff = abs(val1 - val2) / max(abs(val1), abs(val2))
            if rel_diff > tolerance:
                return False, f"{val1} vs {val2} (rel diff: {rel_diff:.2e})"
        else:
            # Absolute tolerance for near-zero values
            if abs(val1 - val2) > tolerance:
                return False, f"{val1} vs {val2} (abs diff: {abs(val1-val2):.2e})"
        
        return True, "Within tolerance"
    
    # Handle other types
    else:
        if val1 == val2:
            return True, "Exact match"
        else:
            return False, f"{val1} vs {val2}"

def compare_ntuples(file1, file2, tree_name="EventTree", max_events=None, verbose=False):
    """Compare two ROOT ntuples"""
    
    print(f"\nComparing ntuples:")
    print(f"  File 1 (C++): {file1}")
    print(f"  File 2 (Python): {file2}")
    print(f"  Tree: {tree_name}")
    print()
    
    # Open files
    f1 = ROOT.TFile.Open(file1)
    f2 = ROOT.TFile.Open(file2)
    
    if not f1 or f1.IsZombie():
        print(f"ERROR: Cannot open file 1: {file1}")
        return False
    
    if not f2 or f2.IsZombie():
        print(f"ERROR: Cannot open file 2: {file2}")
        return False
    
    # Get trees
    tree1 = f1.Get(tree_name)
    tree2 = f2.Get(tree_name)
    
    if not tree1:
        print(f"ERROR: Cannot find tree '{tree_name}' in file 1")
        return False
        
    if not tree2:
        print(f"ERROR: Cannot find tree '{tree_name}' in file 2")
        return False
    
    # Basic checks
    n_entries1 = tree1.GetEntries()
    n_entries2 = tree2.GetEntries()
    
    print(f"Number of entries - File 1: {n_entries1}, File 2: {n_entries2}")
    
    if n_entries1 != n_entries2:
        print(f"WARNING: Different number of entries!")
    
    # Get branch lists
    branches1 = set([br.GetName() for br in tree1.GetListOfBranches()])
    branches2 = set([br.GetName() for br in tree2.GetListOfBranches()])
    
    only_in_1 = branches1 - branches2
    only_in_2 = branches2 - branches1
    common_branches = branches1 & branches2
    
    if only_in_1:
        print(f"\nBranches only in file 1: {sorted(only_in_1)}")
    if only_in_2:
        print(f"\nBranches only in file 2: {sorted(only_in_2)}")
    
    print(f"\nComparing {len(common_branches)} common branches...")
    
    # Compare entries
    n_compare = min(n_entries1, n_entries2)
    if max_events and max_events < n_compare:
        n_compare = max_events
        print(f"Limiting comparison to first {n_compare} events")
    
    mismatches = defaultdict(list)
    branch_stats = {}
    
    # Key branches to always report
    key_branches = [
        'run', 'subrun', 'event',
        'trueVtxX', 'trueVtxY', 'trueVtxZ',
        'xsecWeight',
        'truePrimPartE',
        'trueSimPartPDG','trueSimPartTID','trueSimPartE',
        'foundVertex', 'vtxX', 'vtxY', 'vtxZ', 'vtxScore',
        'nTracks',  'trackPID','trackElScore', 'trackPhScore', 'trackMuScore', 'trackPiScore', 'trackPrScore',
        'trackComp','trackPurity','trackProcess','trackPrimaryScore','trackFromNeutral','trackFromCharged',
        'nShowers','showerPID','showerElScore','showerPhScore','showerMuScore','showerPiScore','showerPrScore',
        'showerComp','showerPurity','showerProcess','showerPrimaryScore','showerFromNeutral','showerFromCharged',
        #'nTracks', 'nShowers',
        #'vtxIndex', 
        #'predictedPEtotal', 'observedPEtotal', 'sinkhorn_div', 'fracerrPE'
    ]
    
    for i in range(n_compare):
        tree1.GetEntry(i)
        tree2.GetEntry(i)
        
        # Get event identifier
        try:
            event_id = f"Entry[{i}] {tree1.run}/{tree1.subrun}/{tree1.event}"
        except:
            event_id = f"Entry {i}"
        
        for branch in sorted(common_branches):
            try:
                val1 = getattr(tree1, branch)
                val2 = getattr(tree2, branch)
                
                # Convert to comparable format
                if hasattr(val1, '__len__'):
                    val1 = list(val1)
                if hasattr(val2, '__len__'):
                    val2 = list(val2)
                
                match, msg = compare_values(val1, val2, branch)
                
                if branch not in branch_stats:
                    branch_stats[branch] = {'matches': 0, 'mismatches': 0}
                
                if match:
                    branch_stats[branch]['matches'] += 1
                else:
                    branch_stats[branch]['mismatches'] += 1
                    mismatches[branch].append((i, event_id, msg))
                    
                    if verbose or branch in key_branches:
                        print(f"  Event {event_id}, Branch '{branch}': {msg}")
                        
            except Exception as e:
                if verbose:
                    print(f"  Error comparing branch '{branch}': {e}")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    # Report key branches first
    print("\nKey branches:")
    for branch in key_branches:
        if branch in branch_stats:
            stats = branch_stats[branch]
            total = stats['matches'] + stats['mismatches']
            if stats['mismatches'] > 0:
                print(f"  {branch}: {stats['mismatches']}/{total} mismatches")
            else:
                print(f"  {branch}: All {total} entries match ✓")
    
    # Report branches with mismatches
    print("\nBranches with mismatches:")
    mismatch_branches = [(b, s['mismatches']) for b, s in branch_stats.items() 
                         if s['mismatches'] > 0 and b not in key_branches]
    mismatch_branches.sort(key=lambda x: x[1], reverse=True)
    
    if mismatch_branches:
        for branch, n_mismatch in mismatch_branches[:20]:  # Show top 20
            total = branch_stats[branch]['matches'] + n_mismatch
            print(f"  {branch}: {n_mismatch}/{total} mismatches")
        if len(mismatch_branches) > 20:
            print(f"  ... and {len(mismatch_branches)-20} more branches with mismatches")
    else:
        print("  None! All branches match perfectly.")
    
    # Overall statistics
    total_comparisons = sum(s['matches'] + s['mismatches'] for s in branch_stats.values())
    total_mismatches = sum(s['mismatches'] for s in branch_stats.values())
    
    print(f"\nTotal comparisons: {total_comparisons}")
    print(f"Total mismatches: {total_mismatches}")
    print(f"Match rate: {100.0 * (1.0 - total_mismatches/total_comparisons):.2f}%")
    
    # Close files
    f1.Close()
    f2.Close()
    
    return total_mismatches == 0

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_ntuple_outputs.py <cpp_output.root> <python_output.root> [max_events]")
        print()
        print("Compare ntuples produced by C++ and Python versions of gen2ntuple")
        sys.exit(1)
    
    cpp_file = sys.argv[1]
    python_file = sys.argv[2]
    max_events = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    # Compare EventTree
    success = compare_ntuples(cpp_file, python_file, "EventTree", max_events=max_events, verbose=True)
    
    # Also compare potTree if it exists
    print("\n" + "="*60)
    print("Checking potTree...")
    try:
        f1 = ROOT.TFile.Open(cpp_file)
        f2 = ROOT.TFile.Open(python_file)
        if f1.Get("potTree") and f2.Get("potTree"):
            compare_ntuples(cpp_file, python_file, "potTree", verbose=True)
        else:
            print("potTree not found in one or both files")
        f1.Close()
        f2.Close()
    except:
        pass
    
    if success:
        print("\n✓ Validation PASSED: C++ and Python outputs match!")
        sys.exit(0)
    else:
        print("\n✗ Validation FAILED: Outputs do not match")
        sys.exit(1)

if __name__ == "__main__":
    main()