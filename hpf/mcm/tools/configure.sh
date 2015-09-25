#!/bin/bash

# usage: ./configure.sh DB_DIRECTORY
# DB_DIRECTORY  -   indicates where script should set Rosetta to look for DBs and configs
# Decompresses MCM and Rosetta DBs to DB_DIRECTORY as well

ZIPLIST="rosetta++.tbz"
echo "Unpacking rosetta code..."
for f in $ZIPLIST
do
  if [ -f $f -a -r $f ]; then
   #echo "F: $f"
   #echo "Filename: ${f##*/}"
   BASENAME=${f%.*}
   if [ -d "$BASENAME" ]; then
     echo "  $BASENAME exists"
   else
     tar -xjvf "$f"
   fi
  else
   echo "Error: Cannot read $f"
  fi
done

# Unpacking databases into DB_DIRECTORY
echo "Decompressing databases to $1 ..."
DBLIST="mammothDbCEMS.tbz rosetta_database.tbz"
if [! -d $1]; then
  echo "$1 is not a valid directory. Making.."
  mkdir $1
fi
for d in $DBLIST
do
  tar xjf $d -C $1
done

# This will configure the db path files with the first argument
OLD="DB_HOME"
NEW=$1
PATHLIST="./paths.txt.tmpl ./list.mammoth.tmpl"
echo "Creating path files..."
for f in $PATHLIST
do
  if [ -f $f -a -r $f ]; then
   #echo "F: $f"
   #echo "Filename: ${f##*/}"
   BASENAME=${f%.*}
   if [ -f "$BASENAME" ]; then
     echo "  $BASENAME exists"
   else
     sed "s|$OLD|$NEW|g" "$f" > "$BASENAME"
   fi
  else
   echo "Error: Cannot read $f"
  fi
done

echo "Compiling binaries..."
if [ -f "robetta_cluster_mod" ]; then
  echo "  rosetta_cluster_mod exists..."
else
  echo "gcc robetta_cluster_mod_directory/pbCluster/robetta_cluster_mod.c -o robetta_cluster_mod -lm"
  gcc robetta_cluster_mod_directory/pbCluster/robetta_cluster_mod.c -o robetta_cluster_mod -lm
fi

if [ -f "rosetta" ]; then
  echo "  rosetta exists..."
else
  cd rosetta++ && make gcc
  cd ../
  ln -s rosetta++/rosetta.gcc rosetta
  echo "  !!YOU BETTER MAKE SURE ROSETTA++ IS LINKED CORRECTLY!!"
fi
