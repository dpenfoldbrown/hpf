#!/bin/bash
for ((a=0; a <=25; a++))
do
echo "Replicating compute-0-$a"
ssh compute-0-$a ~/fragment_picker/scripts/replicate.sh
done
