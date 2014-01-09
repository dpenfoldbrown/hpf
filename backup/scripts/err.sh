#!/bin/bash

LOG=$HOME/backup/log
echo `date` > $LOG/err_backup.log
rsync -av -e ssh --exclude '*.md5' dpb@err.bio.nyu.edu:/var/archives/ /archive/rb133/backup/err/ >> $LOG/err_backup.log 2>&1
mailx -s "$? ERR Archive" dpb3@nyu.edu < $LOG/err_backup.log

