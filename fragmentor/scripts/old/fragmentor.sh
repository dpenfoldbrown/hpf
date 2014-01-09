#!/bin/sh
export SGE_TASK_ID=$SGE_TASK_ID
hostname
echo $SGE_TASK_ID
python ./src/frag_array.py ~/frags $1 ~/.local/share
