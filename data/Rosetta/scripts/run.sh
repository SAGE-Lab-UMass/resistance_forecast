#!/bin/sh
# Provide the last residue number
end=11
# Provide the file name (filename="relaxed_200_0.pdb")
filename="output_300_3.pdb"
# You can change this value if want to use more cpus
increment=10
# Change this value if you do not want to start from residue 1
start=1
int_end=$increment

while [ $start -le $end ]; do
  python ddG_calculation.py $start $int_end $filename $end
  echo $start $int_end $filename $end
  start=$((start + increment))
  int_end=$((int_end + increment))
done
