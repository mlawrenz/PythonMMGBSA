#  !/bin/bash
#BSUB -J test-md           # job name
#BSUB -W 72:00                # wall-clock time (hrs:mins)
#BSUB -n 8                  # number of tasks in job
#BSUB -q quickie              # queue
#BSUB -R span[hosts=1]
#BSUB -B
#BSUB -eo test-md.err       # error file name in which %J is replaced by the job ID
#BSUB -oo test-md.out      # output file name in which %J is replaced by the job ID
 
ligfile=mol2-files/2683187.mol2
protfile=pdb-files/pde10a_2683187_protonly.pdb

# implicit minimization only
RunMMGBSA.py -time -jname test1 -prot \
$protfile -im -prad 5.0 -netc 0 -mol2 \
$ligfile  >& test1.log

# explicit minimization only
RunMMGBSA.py -time -jname test2 -prot \
$protfile -prad 5.0 -netc 0 -mol2 \
$ligfile  >& test2.log


# explicit MD, TIP3P water
# default MD step # is 100000 (200 ps)
RunMMGBSA.py -time -jname test3 -prot \
$protfile -prad 5.0 -netc 0 -mol2 $ligfile -md -nproc 8 -mdsteps 10000 >& test3.log


# implicit MD (not actually faster in all cases)
# default MD step # is 100000 (200 ps)
RunMMGBSA.py -time -jname test4 -prot \
$protfile -prad 5.0 -netc 0 -mol2 $ligfile -md -nproc 8 -mdsteps 10000 -im >& test4.log


# explicit GPU MD: run via gpu_calc queue
# assumes 4 possible gpus on ussf-papp-gpu01, will pick an available
# one
RunMMGBSA.py -time -jname test5 -prot \
$protfile -prad 5.0 -netc 0 -mol2 $ligfile -md -nproc 8 -mdsteps 10000 -gpu >& test6.log

