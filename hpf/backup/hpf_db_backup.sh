#!/bin/bash

# Simple mysqldump-based backup for databases

BACKUP_DIR="/data/backup/hpf_db/"
DB_USER="dpb"
DB_PWD="dpb_nyu"
DATABASES="hpf"

cd $BACKUP_DIR

for D in $DATABASES
do
    echo "+ backing up MySQL DB $D to $BACKUP_DIR"
    # Note: must not have space after -p, or mysql prompts for pwd
    mysqldump -u $DB_USER -p$DB_PWD $D | bzip2 -c > $D.sql.bz2
done
echo "Complete."
