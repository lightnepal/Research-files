#!/bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sharedapps/LS/binary/matlab_R2014a/runtime/glnxa64:/sharedapps/LS/binary/matlab_R2014a/bin/glnxa64:/sharedapps/LS/binary/matlab_R2014a/sys/os/glnxa64/
./mint7 1 "$SLURM_ARRAY_TASK_ID"




