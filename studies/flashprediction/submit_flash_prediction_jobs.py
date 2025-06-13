#!/usr/bin/env python3
"""
Script to submit flash prediction jobs to the grid/cluster

This script:
1. Takes a list of dlmerged and reco file pairs
2. Splits them into manageable job chunks
3. Creates submission scripts for each job
4. Submits to SLURM or other batch system
"""

import os
import sys
import argparse
import glob
from pathlib import Path

def create_slurm_script(job_name, dlmerged_file, reco_file, output_file, 
                       num_entries=-1, start_entry=0, log_dir="logs"):
    """Create a SLURM submission script for a single job"""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8GB
#SBATCH --partition=batch
#SBATCH --account=tufts

# Setup environment
echo "Setting up environment..."
cd $WORK_DIR
source /cluster/tufts/wongjiradlabnu/containers/setenv_py3.sh
source /cluster/tufts/wongjiradlabnu/twongj01/ubdl/configure.sh

# Navigate to flash prediction directory
cd /cluster/tufts/wongjiradlabnu/twongj01/gen2ntuple/studies/flashprediction

# Run the executable
echo "Running flash prediction..."
./build/flashprediction/calculate_flash_predictions \\
    -d {dlmerged_file} \\
    -r {reco_file} \\
    -o {output_file} \\
    -n {num_entries} \\
    -s {start_entry} \\
    -v

echo "Job completed!"
"""
    
    return script_content

def find_matching_files(dlmerged_pattern, reco_pattern):
    """Find matching dlmerged and reco file pairs"""
    
    dlmerged_files = sorted(glob.glob(dlmerged_pattern))
    reco_files = sorted(glob.glob(reco_pattern))
    
    # Match files based on run/subrun/event identifiers
    matched_pairs = []
    
    for dlmerged in dlmerged_files:
        # Extract identifier from filename (customize based on your naming convention)
        base_name = os.path.basename(dlmerged)
        # Example: look for matching reco file
        for reco in reco_files:
            if base_name.replace("dlmerged", "reco") in reco:
                matched_pairs.append((dlmerged, reco))
                break
    
    return matched_pairs

def main():
    parser = argparse.ArgumentParser(description="Submit flash prediction jobs to cluster")
    
    parser.add_argument("-d", "--dlmerged-pattern", required=True,
                        help="Glob pattern for dlmerged files (e.g., '/path/to/dlmerged*.root')")
    parser.add_argument("-r", "--reco-pattern", required=True,
                        help="Glob pattern for reco files (e.g., '/path/to/reco*.root')")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for flash prediction files")
    parser.add_argument("-j", "--jobs-per-file", type=int, default=1,
                        help="Number of jobs per file pair (for splitting large files)")
    parser.add_argument("-e", "--entries-per-job", type=int, default=-1,
                        help="Number of entries per job (default: all)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Create scripts but don't submit")
    parser.add_argument("--job-dir", default="jobs",
                        help="Directory to store job scripts")
    parser.add_argument("--log-dir", default="logs",
                        help="Directory for job logs")
    
    args = parser.parse_args()
    
    # Create necessary directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.job_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)
    
    # Find matching file pairs
    print("Finding matching file pairs...")
    file_pairs = find_matching_files(args.dlmerged_pattern, args.reco_pattern)
    
    if not file_pairs:
        print("No matching file pairs found!")
        sys.exit(1)
    
    print(f"Found {len(file_pairs)} file pairs")
    
    # Create job scripts
    job_count = 0
    
    for dlmerged_file, reco_file in file_pairs:
        base_name = Path(dlmerged_file).stem
        
        if args.jobs_per_file == 1:
            # Single job per file
            job_name = f"flash_pred_{base_name}"
            output_file = os.path.join(args.output_dir, f"flash_predictions_{base_name}.root")
            
            script_content = create_slurm_script(
                job_name, dlmerged_file, reco_file, output_file,
                num_entries=args.entries_per_job, start_entry=0,
                log_dir=args.log_dir
            )
            
            script_path = os.path.join(args.job_dir, f"{job_name}.sh")
            with open(script_path, 'w') as f:
                f.write(script_content)
            
            os.chmod(script_path, 0o755)
            
            if not args.dry_run:
                os.system(f"sbatch {script_path}")
                print(f"Submitted: {job_name}")
            else:
                print(f"Created script: {script_path}")
            
            job_count += 1
            
        else:
            # Multiple jobs per file (for large files)
            # First, we'd need to determine total entries - simplified here
            entries_per_job = args.entries_per_job if args.entries_per_job > 0 else 1000
            
            for job_idx in range(args.jobs_per_file):
                start_entry = job_idx * entries_per_job
                job_name = f"flash_pred_{base_name}_part{job_idx}"
                output_file = os.path.join(args.output_dir, 
                                         f"flash_predictions_{base_name}_part{job_idx}.root")
                
                script_content = create_slurm_script(
                    job_name, dlmerged_file, reco_file, output_file,
                    num_entries=entries_per_job, start_entry=start_entry,
                    log_dir=args.log_dir
                )
                
                script_path = os.path.join(args.job_dir, f"{job_name}.sh")
                with open(script_path, 'w') as f:
                    f.write(script_content)
                
                os.chmod(script_path, 0o755)
                
                if not args.dry_run:
                    os.system(f"sbatch {script_path}")
                    print(f"Submitted: {job_name}")
                else:
                    print(f"Created script: {script_path}")
                
                job_count += 1
    
    print(f"\nTotal jobs: {job_count}")
    if args.dry_run:
        print("Dry run complete - no jobs submitted")
        print(f"To submit, run without --dry-run")

if __name__ == "__main__":
    main()