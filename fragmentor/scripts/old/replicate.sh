#!/bin/bash
TMP='/state/partition1/bonneaulab/fragmentor'
mkdir -p $TMP
echo "  Removing $TMP/*"
rm -rf $TMP/*
echo "  Untar"
tar -C $TMP -xzvf ~/fragment_picker/data.tgz
