#! /bin/bash
#
# Bash script to submit LAMMPS simulation
# on local PC. 
#
# Input/Output
NPROC=1
INPUT_FOLDER="${PWD}"
INPUT_FILE="NVT_class2xe.script"
OUTPUT_FILE="out.out"

# LAMMPS install/build location
LAMMPS="/mnt/c/Users/jdkem/Desktop/GitHub/LAMMPS/lammps_class2xe"

# Run using the binary from the build directory
export OMP_NUM_THREADS=1
mpirun -n   "${NPROC}" "${LAMMPS}/build/lmp" \
       -log "${INPUT_FOLDER}/${OUTPUT_FILE}" \
       -in  "${INPUT_FOLDER}/${INPUT_FILE}"