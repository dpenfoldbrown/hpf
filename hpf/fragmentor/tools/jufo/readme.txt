JUFO secondary structure readme 04/26/2002
contact jens@jens-meiler.de for questions

The directory does contain the executable decoys_protein_sec.exe, 
that performs the actual secondary structure preditcion. Since
a scoring matrix from PSIBLAST is necessary for all predictions,
the shell scripts will perform a PSIBLAST run and call JUFO
afterwards. Three scripts are available:

* "run_jufo_1D" needs only the "*.fasta" file and predicts secondary 
  structure form sequence only. The output is written to "*.jufo_1D_ss". 
  It will report TWO warnings, since no PDB file and no ROSETTA silent 
  mode file could be found => dont care.

* "run_jufo_3Dpdb" needs the "*.fasta" file and a "*.pdb" file and predicts 
  secondary structure using the sequence and low resolution structural 
  information from the PDB file. The output is written to "*.jufo_3D_ss".
  It will report ONE warning, since no ROSETTA silent mode file could
  be found => dont care. The protocol will also perform the sequence only 
  prediction and write the data to "*.jufo_1D_ss".

* "run_jufo_3Ddecoys" needs the "*.fasta" file and a "*.out" file and predicts
  secondary structure using the sequence and low resolution structural
  information from all decoys in the OUT file. The output is written to
  "*.decoys.jufo_3D_ss". It will report ONE warnings, since no PDB file 
  could be found => dont care. The protocol will also perform the sequence only
  prediction and write the data to "*.jufo_1D_ss".
