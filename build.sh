#!/bin/bash
Green='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color
mpicxx "$1.cpp" -o "$1" -lm
printf "${Green}Finished Compiling, Started running${NC}\n"
mpirun -np "$2" "./$1" 