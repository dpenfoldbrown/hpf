#!/bin/bash

LOG=$HOME/backup/log
echo `date` > $LOG/hpf_backup.log
rsync -av -e ssh dpb@markula.bio.nyu.edu:/home/wcgrid/incoming/ /archive/rb133/hpf2/results >> $LOG/hpf_backup.log 2>&1
mailx -s "$? HPF Archive" dpb3@nyu.edu < $LOG/hpf_backup.log

