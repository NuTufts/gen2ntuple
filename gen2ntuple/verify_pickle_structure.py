#!/usr/bin/env python3
"""
Verify the structure of the weight pickle files to ensure they match
what the C++ code expects: dict[run][subrun][event] = weight
"""

import sys
import pickle

def verify_pickle_structure(pkl_file):
    """Load and verify the structure of a weight pickle file"""
    
    print(f"Loading pickle file: {pkl_file}")
    
    with open(pkl_file, 'rb') as f:
        weights = pickle.load(f)
    
    # Check that it's a dictionary
    if not isinstance(weights, dict):
        print(f"ERROR: Root object is not a dict, it's a {type(weights)}")
        return False
    
    print(f"Total number of runs: {len(weights)}")
    
    # Sample a few entries to verify structure
    sample_count = 0
    total_events = 0
    
    for run, run_data in weights.items():
        if not isinstance(run_data, dict):
            print(f"ERROR: Data for run {run} is not a dict")
            return False
            
        for subrun, subrun_data in run_data.items():
            if not isinstance(subrun_data, dict):
                print(f"ERROR: Data for run {run} subrun {subrun} is not a dict")
                return False
                
            for event, weight in subrun_data.items():
                total_events += 1
                
                # Print first few examples
                if sample_count < 5:
                    print(f"  Example: run={run}, subrun={subrun}, event={event}, weight={weight}")
                    sample_count += 1
                
                # Verify weight is a number
                if not isinstance(weight, (int, float)):
                    print(f"ERROR: Weight for run {run} subrun {subrun} event {event} is not a number: {type(weight)}")
                    return False
    
    print(f"Total number of events with weights: {total_events}")
    print("Pickle file structure verified successfully!")
    
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <weight_file.pkl>")
        sys.exit(1)
    
    pkl_file = sys.argv[1]
    
    if verify_pickle_structure(pkl_file):
        sys.exit(0)
    else:
        sys.exit(1)