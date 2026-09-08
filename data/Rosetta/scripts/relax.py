"""
Protein Structure Relaxation Script using PyRosetta

This script performs energy minimization (FastRelax) on a protein structure.
It can either:
1. Determine an optimal number of iterations (`max_iter`) for relaxation, or
2. Perform relaxation a fixed number of times with a user-specified `max_iter`.

Usage:
    python relax.py input_structure.pdb <max_iter>

Arguments:
    input_structure.pdb : Path to the input PDB file
    max_iter            : Maximum number of iterations for relaxation
                          (set to 0 to automatically find optimal max_iter)
"""

import sys
import argparse
from pyrosetta import *
from rosetta.protocols.relax import FastRelax

# Initialize PyRosetta
init()

# Create the full-atom scoring function
scorefxn = get_fa_scorefxn()

# Create argument parser
parser = argparse.ArgumentParser(description="Protein Structure Relaxation using PyRosetta (FastRelax).")

# Add arguments
parser.add_argument("--input_pdb", type=str, help="Path to the input PDB file to be relaxed.")
parser.add_argument("--max_iter", type=int, help="Maximum number of iterations for relaxation. Set to 0 to automatically find optimal iterations.")

# Parse arguments
args = parser.parse_args()
input_pdb = args.input_pdb
max_iter = args.max_iter

# Initialize tracking of lowest energy structure
pose = pose_from_pdb(input_pdb)
lowest_energy = scorefxn.score(pose)
lowest_energy_pdb = input_pdb

# If max_iter is 0, automatically determine optimal iterations
if max_iter == 0:
    print("Finding optimal max_iter...")
    for i in [100, 200, 300, 400, 500]:
        pose = pose_from_pdb(input_pdb)
        relax = FastRelax()
        relax.set_scorefxn(scorefxn)
        relax.max_iter(i)
        relax.apply(pose)
        output_file = f'output_{i}.pdb'
        pose.dump_pdb(output_file)

        energy = scorefxn.score(pose)
        if energy <= lowest_energy:
            lowest_energy = energy
            lowest_energy_pdb = output_file
            max_iter = i

# Perform relaxation 4 times using determined or provided max_iter
print(f"Performing relaxation 4 times with max_iter = {max_iter}...")
for iteration in range(4):
    pose = pose_from_pdb(input_pdb)
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.max_iter(max_iter)
    relax.apply(pose)

    output_file = f'output_{max_iter}_{iteration}.pdb'
    pose.dump_pdb(output_file)

    energy = scorefxn.score(pose)
    if energy <= lowest_energy:
        lowest_energy = energy
        lowest_energy_pdb = output_file

# Summary
print(f"\nLowest energy pose: {lowest_energy_pdb}")
print(f"Energy: {lowest_energy}")
print(f"Max iterations used: {max_iter}")
