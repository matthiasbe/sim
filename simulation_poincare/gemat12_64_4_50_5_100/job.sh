#!/bin/bash
#@ class            = clallmds+
#@ job_name         = MIS_BG
#@ total_tasks      = 64
#@ node             = 4
#@ node_usage       = not_shared
#@ wall_clock_limit = 00:05:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ queue
#
module load intel openmpi
module load mkl/11.0
mpirun ./build/bench -i 50 -e 5 -m 100 -f test/matrices/gemat12.mtx
