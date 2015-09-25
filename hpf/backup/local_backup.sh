#!/bin/bash

# A standard backup script for tar'ing directories and SQL
# Config:

BACKUP_USER="hpf-backup"
BACKUP_HOST="handbanana.bio.nyu.edu"
BACKUP_KEY="/home/dpb/.ssh/remote_backup.rsa"
REMOTE_DIR="/data/backup/markula/"

BACKUP_DIR="/data/backup"
cd $BACKUP_DIR

echo "Backup started: " `date`


## Directory section - tars directories to back up and puts them in
# directory BACKUP_DIR. 

DIR="/var/webapp"
NAME="markula-var-webapp.tar.gz"

echo "Backing up $DIR to $BACKUP_DIR"
tar czf $NAME $DIR



## Database section - sqldumps databases, tars the dump, and puts
# them in BACKUP_DIR. Add space-delim'ed entries to DATABASES
# for more

DATABASES="grid"

for D in $DATABASES
do
        echo "Backing up MySQL DB $D to $BACKUP_DIR"
        mysqldump --databases $D | bzip2 -c > $D.sql.bz2
done

## Or all DBs
#echo "Backing up all MySQL DBs to $BACKUP_DIR"
#mysqldump --all-databases | bzip2 -c > all-databases.sql.bz2



## Rsync section - rsyncs all backups to backup host
#
echo "Sync'ing backups to $BACKUP_HOST:$REMOTE_DIR"
rsync -av -e "ssh -i $BACKUP_KEY" --exclude 'local_backup.*' $BACKUP_DIR $BACKUP_USER@$BACKUP_HOST:$REMOTE_DIR



echo "Backup complete: " `date`

