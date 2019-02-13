#!/bin/bash
#@ class            = clallmds+
#@ job_name         = MIS_BG
#@ total_tasks      = 1
#@ node             = 1
#@ node_usage       = not_shared
#@ wall_clock_limit = 00:05:00
#@ output           = $(job_name).log
#@ error            = $(job_name).err
#@ job_type         = mpich
#@ environment      = COPY_ALL
#@ queue
#
module load intel openmpi
module load mkl/11.0
export OMP_NUM_THREADS=16
mpirun -x OMP_NUM_THREADS -bind-to-socket -report-bindings ./build/bench -i 50 -e 5 -m 10 -f test/matrices/gemat12.mtx
