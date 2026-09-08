# Structure Relaxation and ddG Calculation Pipeline

This pipeline performs structure relaxation and calculates ddG values using PyRosetta for all possible variants at each residue position. Follow the steps below to run the necessary scripts.

## Steps

1. **Relax the Structure:**
   - Use the `relax.py` script to relax the initial PDB structure with REF2015 score function.
   - Provide the initial PDB file as input.
   - Choose an optimized `max_iter` value for the relaxation process.
   - Run `relax.py` with the optimized `max_iter` 4 times to obtain the best-optimized structure.
   ```
   python relax.py initial_structure.pdb 0
   ```

2. **Select the Lowest Energy Pose:**
   - After performing relaxation, select the lowest energy pose from the output structures.
   - This pose will be used as input for ddG calculations.

3. **Run ddG Calculations:**
   - The `run.sh` script will automatically call `ddG_calculation.py` for all possible variants at each position.
   - Provide the lowest energy PDB structure and the total number of residues as input.

## Notes
- **Residue Index Mismatch:**
  - The residue index in the original PDB file may differ from Rosetta numbering.
  - If needed, use `convert_ros2pdb.py` to convert between Rosetta and PDB residue numbering.
      ```
      python convert_ros2pdb.py "variant.csv" "new_file.csv" 0 1
      ```

## Example Usage
```
# Relax structure with 0 max iterations to find optimal max_iter
python relax.py initial_structure.pdb 0

# Run ddG calculations
./run.sh
```