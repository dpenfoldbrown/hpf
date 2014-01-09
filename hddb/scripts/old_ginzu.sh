#!/bin/bash
echo $CWD
source $HOME/.bashrc
./washington.sh
./ginzu.py $1
exit $?