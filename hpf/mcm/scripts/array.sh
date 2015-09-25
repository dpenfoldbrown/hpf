#!/bin/bash
#PBS -l nodes="1:ppn=2"
#PBS -l walltime=20:00:00
#PBS -q ser2
#PBS -S /bin/bash
#PBS -V 
#PBS -j oe 
#PBS -m be 
#PBS -o /scratch/pcw216/tmp/mcm/logs

#echo "HOSTNAME= $HOSTNAME"
#echo "SHELL=$SHELL"
#echo "PATHS= $PATH"
#echo "PWD= $PWD"
#echo "PYTHONPATH= $PYTHONPATH"
TMP="/scratch/pcw216/tmp"
SRC=/home/pcw216/local/projects/mcm
mkdir -p "$TMP"
# LOCK="$TMP/mcmlock"
# echo "LOCK defined ${LOCK}"
# while [ -e ${LOCK} ]; do
# 	MYPID=`head -n 1 "${LOCK}"`
# 	TEST_RUNNING=`ps -p ${MYPID} | grep ${MYPID}`
# 	if [ -z "${TEST_RUNNING}" ]; then
# 		rm -f ${LOCK}
# 	else
# 		sleep 15
# 	fi
# done
# echo "Touching ${LOCK}"
# echo $$ > ${LOCK}

DIRECTORY="$TMP/mcm"
if [ -d ${DIRECTORY} ]; then
	echo "MCM exists"
else
	echo "MCM doesn't exist"
	mkdir -p ${DIRECTORY}
fi

DIRECTORY="$TMP/mcm/db"
if [ -d ${DIRECTORY} ]; then
	echo "DB exists"
else
	echo "DB doesn't exist"
	exit 1
#	./scripts/replicate.sh $TMP $SRC
fi
#echo "Removing ${LOCK}"
#rm -f ${LOCK}

# 2,3 are for any alternative options
#CODE="lm"
#TASK=1
RESULTS=$TMP/mcm/results/$CODE
OUTFILE=$TMP/mcm/input/$CODE/$CODE\_hpf2_results/$CODE%s.result.gz
PSIPRED=$TMP/mcm/input/$CODE/$CODE%s/$CODE%s.psipred

cd $SRC
echo "python ./src/mcm_array.py -t $TASK -s 15000 -g 10000 -d $TMP/mcm -r $RESULTS $CODE $OUTFILE $PSIPRED"
python ./src/mcm_array.py -t $TASK -s 15000 -g 10000 -d $TMP/mcm -r $RESULTS $CODE $OUTFILE $PSIPRED
exit $?
