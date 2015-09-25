#!/usr/bin/env python

## import_fasta
## An object script to read and parse a given fasta file (multiple formats
## supported, see below), and to import all protein sequences found in the file
## to a database.
##
## Duncan Penfold-Brown, 2/15/2011
##
## Uses BioPython (Bio.SeqIO) to parse fasta records in the default FastaFile
##  object.
## Database connectivity and ORM managed in hpf.hddb.db - database specifcations
##  are defined therein (including db location, user, pass, db to use, etc).
## Extend the FastaFile class to support advanced/custom parsing of fasta files
##  (eg., create class TIGRFastaFile(FastaFile) with custom parsing support methods)
## 
## Some fasta formats: http://en.wikipedia.org/wiki/FASTA_format

## Imports
import os
import re
import sys
from Bio import SeqIO
from hashlib import sha1
from datetime import datetime
from optparse import OptionParser
from sqlalchemy.exc import IntegrityError 
# hpf.hddb.db for real, hpf.hddb.db_test for fakesies (hpf_dev DB).
from hpf.hddb.db import Session, engine, Experiment, Sequence, Protein, SequenceAc, MouseIDMap


## Globals

usage_str = """Usage: %prog [-f | --format FORMAT] [-d | --source_db STRING] [-p | --protein_type STRING] <fasta_file> <experiment_id>

    <fasta_file>    Fasta file to be parsed into records that will be stored in DB.
    <experiment_id> The ID of the experiment that fasta records read from input
                      file will be linked to (must be an integer and valid ID for
                      pre-created experiment table).
    
    Supported formats (-f tag): NCBI (ncbi standard), TIGR, UNIPROT, MOUSE, REFSEEK_CUSTOM, PHYTOZOME.

    [-p | --protein_type] is simply a value that will be stored in the 'protein_type' field of the hpf.protein DB table

    Note: DB specifications in hpf.hddb.db module
    
    Warning: import_fasta will not create duplicate sequences in the database. If a sequence in the database matches a sequence
    in the input file, new protein and sequenceAc table for the sequence will be linked to the preexisting DB sequence.
    This program will, however, create duplicate (identical except for ID) protein and sequenceAc entries in the DB if duplicate
    sequences exist in the input file (ie, multiple protein entries may have the same values and map to the same sequences).
"""


## FastaFile classes

class FastaFile():
# The base class for representing fasta input files and parsing them into individual protein sequences.
# Uses a basic, unsophisticated parse via Bio.SeqIO, producing SeqRecord objects.
# SUBCLASS parsing methods to refine parsing, and database methods to support adding custom info to
# the database. (Generally, _refine_records() and _push_sequenceAc() are the two methods to override).
# Parameters: 
#   fasta_file (filename string)
#   experiment_id (ID of manually created experiment record in database, string of digits) 
#   source_db (The database from which sequences were taken, string).
#   protein_type (A description of the type of protein sequences represent, string).

    def __init__(self, fasta_file, experiment_id, source_db="na", protein_type="bioinformatics"):
        self.fasta_file = fasta_file
        self.experiment_id = experiment_id
        self.source_db = source_db
        self.protein_type = protein_type
        self.records = None
        self.experiment = None
        
        # Create variable to hold SQLAlchemy db session. (Populated when needed).
        self.session = None
        
        # Create count variables for parse and import statistics.
        self.num_records = 0
        self.num_importedrecords = 0
        self.num_skippedseqs = 0

        print "Creating FastaFile object."


    def parse(self, ):
        print "FastaFile parse."
       
        # Open given fasta file and parse fasta.
        handle = open(self.fasta_file, "rU")
        self.records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        
        # Do extra parsing on basic BioPython SeqRecords (store extra values in record.annotation dict).
        self._refine_records()
        
        # Return list of parsed sequence records
        self.num_records = len(self.records)
        print "Number of records parsed from {0}: {1}".format(self.fasta_file, self.num_records)
        return self.records
        
        
    def _refine_records(self, ):
    # Parse remaining required values out of the information in the SeqRecords returned by SeqIO.parse.
    # Store additional values (accession_id, etc) in the SeqRecord.annotations dict, in [key, value] pairs
    # where key is the name of the field (eg. accession_id, gi_number, source_db, etc).
    # Subclass THIS method in order to support particular additional fields for NCBI, TIGR, etc formats.
    # NOTE: acts directly on instance variable self.records (no passing required).
    # NOTE: introduces some overhead, as we're now iterating over the entire list of records twice. Blah.
        
        print "FastaFile _refine_records."
        
        ac_pattern = "(?P<accession_id>[a-zA-Z]{2}[0-9_.]+)"
    
        for record in self.records:
            
            # Best-effort find the accession id of the parsed sequence. If no match, set accession to "".
            ac_found = re.search(ac_pattern, record.id)
            if ac_found:
                record.annotations['accession_id'] = ac_found.group('accession_id')
            else:
                record.annotations['accession_id'] = ""
                
                             
    def db_import(self, ):
    # Pushes all sequences in self.records to the database (spec'ed in hpf.hddb.db).
    # Puts entries in the following tables: Sequence, Protein, SequenceAc
    # Protein links: Sequence (id), Experiment (id) 
    # SequenceAc links: Sequence (id), Protein (id), Experiment (taxonomy_id)
    # NOTE: If pushing a sequence to the DB fails, program will exit with the input file partially
    # uploaded to the DB. Must manually split file and import remaining sequences after error is solved.
        print "FastaFile db_import."
        
        # Initialize DB if no session exists.
        if self.session == None:
            self._db_setup()
        
        # For each record from given fasta file, create ORM objects for all tables and push them to DB.
        for record in self.records:
            
            # Create and push Sequence, Protein, and SequenceAc objects via support methods.
            try:
                seq = self._push_sequence(record)
                prot = self._push_protein(record, seq)
                seqAc = self._push_sequenceAc(record, seq, prot)
            except:
                print "Last record attempted: "
                print record
                print "(This sequence has NOT been added to the DB. All prior sequences in input file have been added to the DB."
                print "Please manually split the file and re-run script with the file containing this sequence and all following sequences."
                print "You may need to manually remove partial records regarding this sequence from the sequence, protein, or sequenceAc tables.)"
                raise
            
            # After pushing a record iteration to the DB, clear the session to avoid object conflicts.
            self.session.expunge_all()
            self.num_importedrecords += 1
        
        # Print import stats.
        print "Total number of imported records: {0}".format(self.num_importedrecords)
        print "Number of sequences skipped (already in DB {0}) : {1}".format(engine.url, self.num_skippedseqs)

    
    # The following _push methods are created to make subclassing to support custom
    # DB insertions/object creations easier.

    def _push_sequence(self, record):
    # Creates a Sequence object, pushes it to the DB, refreshes it (to update id value), returns updated object.
        seq_obj = Sequence(sequence=str(record.seq), sha1=sha1(str(record.seq)).hexdigest())
        self.session.add(seq_obj)
        
        # Try to push sequence to DB - if the same sequence is already in DB, set seq object to preexisting
        # sequence and link latter entries to it (don't want to and cannot duplicate sequences in DB).
        try: 
            self.session.flush()
        except IntegrityError:
            print "Sequence {0}, {1} already in database. Setting sequence object to preexisting sequence.".format(seq_obj.sha1, seq_obj.sequence)
            self.session.rollback()
            seq_obj = self.session.query(Sequence).filter_by(sha1=seq_obj.sha1).first()
            self.num_skippedseqs += 1
            print "Set sequence object to: "
            print seq_obj
        except:
            print "Error pushing sequence to sequence table in DB."
            raise
                        
        self.session.refresh(seq_obj)
        return seq_obj
    
    def _push_protein(self, record, sequence):
    # Creates, pushes, and refreshes a Protein object. Takes sequence record and updated sequence object.
    # Links to experiment via 'id' and to sequence via 'id'.
        prot_obj = Protein(experiment_key=self.experiment_id, protein_type=self.protein_type, sequence_key=sequence.id, probability=0, insert_date=datetime.now())
        self.session.add(prot_obj)
        try:
            self.session.flush()
        except:
            print "Error pushing protein to protein table in DB."
            raise
        self.session.refresh(prot_obj)
        return prot_obj

    def _push_sequenceAc(self, record, sequence, protein):
    # Creates, pushes, and refreshes a SequenceAc obj. Takes sequence record, sequence object, and protein object.
    # Links to sequence (via id), protein (via id), and experiment (via taxonomy_id).
    # Note that record.annotations dict is populated in the class method _refine_records.
        seqAc_obj = SequenceAc(sequence_key=sequence.id, protein_key=protein.id, db=self.source_db, ac=record.annotations['accession_id'], ac2="", description=record.description, taxonomy_id=self.experiment.taxonomy_id, insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj

    def _db_setup(self, ):
        print "FastaFile db_setup"
        
        # Create a session via hpf.hddb.db Session (a Session object from sessionmaker, sqlalchemy).
        self.session = Session()
        print "Connected to database: {0}".format(engine.url)
        
        # Query the DB session for the given experiment table (if none, no experiment of that ID).
        self.experiment = self.session.query(Experiment).filter(Experiment.id == self.experiment_id).first()
        if self.experiment == None: 
            raise ValueError("Experiment {0} does not exist in the database {1}.".format(self.experiment_id, engine.url))
        
class NCBIFastaFile(FastaFile):
# Fasta header format (NCBI standard):
# >gi|83316026|ref|XP_731047.1| ribosomal protein L10 [Plasmodium yoelii yoelii str. 17XNL]
# >gi | gi number | database | accession ID | description (locus)
# (http://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml)

    def __init__(self, fasta_file, experiment_id, source_db="na", protein_type="bioinformatics"):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating NCBIFastaFile object."

    def _refine_records(self, ):
    # Custom record population for NCBI fasta headers. Parses sequence GI number, source db, and
    # accession id. Also, the description field of the record is fixed (pared down from whole header).  
        print "NCBIFastaFile _refine_records."
                
        ncbi_pattern = r"gi\|(?P<gi_number>\d+)\|(?P<source_db>[a-zA-Z0-9]{2,5})\|(?P<accession_id>[a-zA-Z0-9_.]+)\|\s*(?P<desc>.*)"

        for record in self.records:
            
            # Best-effort find gi_number, source_db, and accession_id in the NCBI fasta header.
            ncbi_found = re.match(ncbi_pattern, record.description)
            if ncbi_found:
                record.annotations['gi_number'] = int(ncbi_found.group('gi_number'))
                record.annotations['source_db'] = ncbi_found.group('source_db')
                record.annotations['accession_id'] = ncbi_found.group('accession_id')
                
                # Fix description field (cut off ID - capture everything in desc. that is not ID).
                record.description = ncbi_found.group('desc')
            else:
                # Set all extra values to defaults and leave description as parsed by Bio.SeqIO.
                record.annotations['gi_number'] = 0
                record.annotations['source_db'] = self.source_db
                record.annotations['accession_id'] = ""
            

    def _push_sequenceAc(self, record, sequence, protein):
    # Creates, pushes, and refreshes a SequenceAc obj. Overrides superclass method to populate DB
    # with information from custom NCBI record parsing. Fields gi, db, ac, and description are new.
        seqAc_obj = SequenceAc(sequence_key=sequence.id, protein_key=protein.id, \
                               gi=record.annotations['gi_number'], \
                               db=record.annotations['source_db'], \
                               ac=record.annotations['accession_id'], \
                               ac2="", \
                               description=record.description, \
                               taxonomy_id=self.experiment.taxonomy_id, \
                               insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj


class TIGRFastaFile(FastaFile):
# Fasta header format: 
# >TIGR|PY03082 | organism=Plasmodium_yoelii_yoelii_str._17XNL | product=Ribosomal protein L10, putative | location=AABL01000871:16386-17520(+) | length=224

    def __init__(self, fasta_file, experiment_id, source_db="na", protein_type="bioinformatics"):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating TIGRFastaFile object."


    def _refine_records(self, ):
    # Custom record population for NCBI fasta headers. Parses sequence GI number, source db, and
    # accession id. Also, the description field of the record is fixed (pared down from whole header).
        
        print "TIGRFastaFile _refine_records."
        
        tigr_pattern = r"TIGR\s*\|\s*(?P<accession_id>[a-zA-Z]{2}[0-9]+)\s*\|\s*organism=(?P<organism>.*)\s*\|\s*product=(?P<desc>.*)\|\s*location="
        
        for record in self.records:
            
            # Best-effort capture accession_id and revise description field.
            tigr_found = re.match(tigr_pattern, record.description)
            if tigr_found:
                record.annotations['accession_id'] = tigr_found.group('accession_id')
                record.description = "{0} [{1}]".format(tigr_found.group('desc'), tigr_found.group('organism'))
            else:
                # Set accession_id to blank and leave description as parsed by default.
                record.annotations['accession_id'] = ""
        
        
class UniprotFastaFile(FastaFile):
# A general class for uniprot fasta files (MouseFastFile should inherit from this, really).
# Fast header format:
# >sp|A0A183|LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1

    def __init__(self, fasta_file, experiment_id, source_db="uniprot", protein_type="uniprot-bioinf"):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating UniprotFastaFile object"

    def _refine_records(self, ):
    # Custom record population for uniprot-style headers (see comment above).
    # Parses uniprot_id, source_db, and description from the uniprot header.
        print "UniprotFastaFile _refine_records"
        
        uniprot_pattern = r"(?P<source_db>[a-zA-Z]{2})\|(?P<uniprot_id>[a-zA-Z0-9_-]+)\|(?P<desc>.+)"

        for record in self.records:
            # Best-effort capture accession_id (=uniprot_id) and source_db, and revise description from uniprot header.
            uniprot_found = re.match(uniprot_pattern, record.description)
            if uniprot_found:
                # Try to grab gene name. Quick, best effort
                try:
                    desc = record.description.split("|")[2]
                    gene_name = re.search(r"GN=(?P<genename>[a-zA-Z0-9]+)", desc).group("genename")
                except Exception:
                    gene_name = ""

                record.annotations['source_db'] = uniprot_found.group('source_db')
                record.annotations['accession_id'] = uniprot_found.group('uniprot_id')
                record.annotations['gene_name'] = gene_name
                record.description = uniprot_found.group('desc')
            else:
                # Set all extra values to defaults and leave description as parsed by Bio.SeqIO.
                record.annotations['source_db'] = self.source_db
                record.annotations['accession_id'] = ""
                record.annotations['gene_name'] = ""


    def _push_sequenceAc(self, record, sequence, protein):
    # Creates, pushes, and refreshes a SequenceAc obj. Overrides superclass method to populate DB
    # with information from custom uniprot record parsing. Fields db, ac, and description differ from superclass.
        seqAc_obj = SequenceAc(
            sequence_key=sequence.id, 
            protein_key=protein.id, 
            db=record.annotations['source_db'], 
            ac=record.annotations['accession_id'], 
            ac2=record.annotations['gene_name'], 
            description=record.description, 
            taxonomy_id=self.experiment.taxonomy_id, 
            insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj


class RefseekRNAFasta(FastaFile):
# A header from berlin lab using refseek IDs and long description
# Fasta header format:
# >NP_002810.1 polypyrimidine tract-binding protein 1 isoform a [Homo sapiens].

    def __init__(self, fasta_file, experiment_id, source_db="ref", protein_type="Human RNA"):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating RefseekRNAFasta object"

    def _refine_records(self, ):
    # Parses refseek id as accession_id and description
        refrna_pattern = r"(?P<accession>[A-Z]+_[0-9]+\.[0-9]+)\s+(?P<description>.*)"
        for record in self.records:
            refrna_found = re.match(refrna_pattern, record.description)
            if refrna_found:
                record.description = refrna_found.group('description')
                record.annotations['accession_id'] = refrna_found.group('accession')
            else:
                record.annotations['accession_id'] = record.id


    def _push_sequenceAc(self, record, sequence, protein):
    # Creates, pushes and refreshed the SequenceAc ORM object (to sequenceAc table). Overrides superclass method for custom vals.
        seqAc_obj = SequenceAc(sequence_key=sequence.id, protein_key=protein.id, db=self.source_db, ac=record.annotations['accession_id'], ac2="", description=record.description, taxonomy_id=self.experiment.taxonomy_id, insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj


class PhytozomeFastaFile(FastaFile):
# A header for data from Phytozome (online resource)
# Header format:
# >Cre08.g369800.t1.1|PACid:19858278

    def __init__(self, fasta_file, experiment_id, source_db="phytozome", protein_type=""):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating PhytozomeFastaFile object"

    def _refine_records(self, ):
        #Cre08.g382100.t1.1|PACid:19858277
        phytozome_pattern = r"(?P<accession>[a-zA-Z0-9.]+)\|PACid:(?P<accession2>[0-9]+)"
        for record in self.records:
            phytozome_found = re.match(phytozome_pattern, record.description)
            if phytozome_found:
                record.annotations['ac'] = phytozome_found.group('accession')
                record.annotations['ac2'] = phytozome_found.group('accession2')
                record.description = "ac, phytozome ID; ac2, PACid"
            else:
                record.annotations['ac'] = ""
                record.annotations['ac2'] = ""

    def _push_sequence(self, record):
        seq_obj = Sequence(sequence=str(record.seq).rstrip('*'), sha1=sha1(str(record.seq).rstrip('*')).hexdigest())
        self.session.add(seq_obj)
        try:
            self.session.flush()
        except IntegrityError:
            print "Sequence {0}, {1} already in database. Setting sequence object to preexisting sequence.".format(seq_obj.sha1, seq_obj.sequence)
            self.session.rollback()
            seq_obj = self.session.query(Sequence).filter_by(sha1=seq_obj.sha1).first()
            self.num_skippedseqs += 1
            print "Set sequence object to: "
            print seq_obj
        except:
            print "Error pushing sequence to sequence table in DB."
            raise

        self.session.refresh(seq_obj)
        return seq_obj
    
    def _push_sequenceAc(self, record, sequence, protein):
        seqAc_obj = SequenceAc(sequence_key=sequence.id, protein_key=protein.id, db=self.source_db, ac=record.annotations['ac'], ac2=record.annotations['ac2'], description=record.description, taxonomy_id=self.experiment.taxonomy_id, insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj


class MouseFastaFile(FastaFile):
# A uniprot fasta file (sequences from swissprot and trembl) with custom ID mappings (uniprot id to
# MGI ID to mousefunc geneID - retrieved from a mapping file, parsed into mouseID_dict, of the form
# key: uniprot_id, value: [mgi_id, gene_id]).
# Fasta header format:
# >sp|Q9CQV8|1433B_MOUSE 14-3-3 protein beta/alpha OS=Mus musculus GN=Ywhab PE=1 SV=3

    def __init__(self, fasta_file, experiment_id, source_db="uniprot", protein_type="mouse-gi"):
        FastaFile.__init__(self, fasta_file, experiment_id, source_db, protein_type)
        print "Creating MouseFastaFile object."
        
        self.mouseID_dict = self._parse_mousemap("/Users/dpb/Documents/mousefunc/mouse_key.txt")
    
    
    def _parse_mousemap(self, mouseID_file):
    # Open static mouse ID mapping file and parse into mouseID dictionary: uniprot_id -> [mgi_id, gene_id]
        try:
            mouseID_handle = open(mouseID_file, "rU")
        except:
            print "Mouse ID mapping file '{0}' does not exist or is inaccessible.".format(mouseID_file)
            raise
        
        mouse_dict = {}
        mouse_pattern = r"(?P<uniprot_id>[a-zA-Z0-9]{6})\s+(?P<mgi_id>MGI:[0-9]+|None)\s+(?P<gene_id>G[0-9]+|None)"
        for line in mouseID_handle:
            mouse_found = re.match(mouse_pattern, line)
            if mouse_found:
                mouse_dict[mouse_found.group('uniprot_id')] = [mouse_found.group('mgi_id'), mouse_found.group('gene_id')]
            else:
                print "Line '{0}' does not match format line '{1}'".format(line, mouse_pattern)
        
        return mouse_dict
    
    
    def _refine_records(self, ):
    # Custom record population for uniprot sequences combined with special mouse data (see init).
    # Parses uniprot_id, source_db, and description from the uniprot header, and mgi_id (=accession_id2)
    # and gene_id from the dictionary obtained from the mouseID mapping file.
    
        print "MouseFastaFile _refine_records"
        
        uniprot_pattern = r"(?P<source_db>[a-zA-Z]{2})\|(?P<uniprot_id>[a-zA-Z0-9]{6})\|(?P<desc>.+)"
        
        for record in self.records:
            # Best-effort capture accession_id (=uniprot_id) and source_db, and revise description from uniprot header.
            uniprot_found = re.match(uniprot_pattern, record.description)
            if uniprot_found:
                record.annotations['source_db'] = uniprot_found.group('source_db')
                record.annotations['accession_id'] = uniprot_found.group('uniprot_id')
                record.description = uniprot_found.group('desc')
                
                # Try to set mgi_id and gene_id form mouseID_dict based on found uniprot/accession id
                try:
                    record.annotations['mgi_id'] = self.mouseID_dict[record.annotations['accession_id']][0]
                    record.annotations['gene_id'] = self.mouseID_dict[record.annotations['accession_id']][1]
                except KeyError:
                    print "Uniprot/accession id '{0}' does not exist in mouse mapping dictionary.".format(record.annotations['accession_id'])
                    record.annotations['mgi_id'] = "None"
                    record.annotations['gene_id'] = "None"
                    
            else:
                # Set all extra values to defaults and leave description as parsed by Bio.SeqIO.
                record.annotations['source_db'] = self.source_db
                record.annotations['accession_id'] = ""
                record.annotations['mgi_id'] = "None"
                record.annotations['gene_id'] = "None"
                
    def db_import(self, ):
    # Pushes all sequences in self.records to the database (spec'ed in hpf.hddb.db).
    # Puts entries in the following tables: Sequence, Protein, SequenceAc, mouse_geneID
    # Linkery:
    # Same as for superclass, except for additional table: Mouse_GeneID.sequence_key -> Sequence.id
        print "MouseFastaFile db_import."
        
        # Initialize DB if no session exists.
        if self.session == None:
            self._db_setup()
        
        # For each record from given fasta file, create ORM objects for all tables and push them to DB.
        for record in self.records:
            
            # Create and push Sequence, Protein, and SequenceAc objects via support methods.
            try:
                seq = self._push_sequence(record)
                prot = self._push_protein(record, seq)
                seqAc = self._push_sequenceAc(record, seq, prot)
                mgid = self._push_mouseIDmap(record, seq)
            except:
                print "Last record attempted: "
                print record
                print "(This sequence has NOT been added to the DB. All prior sequences in input file have been added to the DB."
                print "Please manually split the file and re-run script with the file containing this sequence and all following sequences."
                print "You may need to manually remove partial records regarding this sequence from the sequence, protein, or sequenceAc tables.)"
                raise
            
            # After pushing a record iteration to the DB, clear the session to avoid object conflicts.
            self.session.expunge_all()
            self.num_importedrecords += 1
        
        # Print import stats.
        print "Total number of imported records: {0}".format(self.num_importedrecords)
        print "Number of sequences skipped (already in DB {0}) : {1}".format(engine.url, self.num_skippedseqs)


    def _push_sequenceAc(self, record, sequence, protein):
    # Creates, pushes, and refreshes a SequenceAc obj. Overrides superclass method to populate DB
    # with information from custom Mouse record parsing. Fields db, ac, ac2 and description are new.
        seqAc_obj = SequenceAc(sequence_key=sequence.id, protein_key=protein.id, db=record.annotations['source_db'], ac=record.annotations['accession_id'], ac2=record.annotations['mgi_id'], description=record.description, taxonomy_id=self.experiment.taxonomy_id, insert_date=datetime.now())
        self.session.add(seqAc_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing protein information to sequenceAc table in DB."
            raise
        self.session.refresh(seqAc_obj)
        return seqAc_obj


    def _push_mouseIDmap(self, record, sequence):
    # Creates, pushes, and refreshes a Mouse_GeneID obj. Mouse_GeneID holds sequence_key to gene_id to mgi_id
    # mappings.
        mgid_obj = MouseIDMap(sequence_key=sequence.id, gene_id=record.annotations['gene_id'], mgi_id=record.annotations['mgi_id'], uniprot_id=record.annotations['accession_id'])
        self.session.add(mgid_obj)
        try:
            self.session.flush()
        except:
            print "Error in pushing Mouse ID mapping information to mouseIDmap table in DB."
            raise
        self.session.refresh(mgid_obj)
        return mgid_obj


## Option Parsing

def set_options(parser):
# Set the command-line options for a passed-in OptionParser object (via optparse)

    parser.add_option("-f", "--format", action="store", dest="format", type="string", metavar="FORMAT",
                      help="The format of the fasta header, used for better parsing. Currently supported: NCBI, UNIPROT, TIGR, MOUSE, REFSEEK_CUSTOM")
    parser.add_option("-d", "--source_db", action="store", dest="source_db", type="string", metavar="DATABASE",
                      default="na",
                      help="The source DB of the sequences imported. EG: amnh, emb, pdb, genbank, ref..")
    parser.add_option("-p", "--protein_type", action="store", dest="protein_type", type="string", metavar="STRING",
                      default="bioinformatics",
                      help="A descriptive type to label proteins imported. EG: bioinformatics, phylogeny, microbiome, amnh..")


def parse_positional_args(args, parser):
# Reads the positional commandline arguments from args list leftover from option parsing.
# Returns a pair: (fasta_file, experiment_key)

    if len(args) != 2:
        print "Error: incorrect usage."
        parser.print_help()
        sys.exit(1)
    
    fasta_file = args[0]
    experiment_id = args[1]
    
    # Check fasta_file legitimacy.
    if not os.path.isfile(fasta_file):
        raise IOError("{0} is not a valid file. Check existence, permissions.".format(fasta_file))
    
    # Check to see if experiment_id is valid (can be treated as an integer).
    try:
        int(experiment_id)
    except:
        parser.print_help()
        raise ValueError("Invalid value for experiment_id. Must be an integer: {0}".format(experiment_id))
    
    return (fasta_file, experiment_id)


## Main

def main():

    # Set and read commandline options and arguments.
    parser = OptionParser(usage=usage_str)
    set_options(parser)
    (options, args) = parser.parse_args()
    (fasta_file, experiment_id) = parse_positional_args(args, parser) 

    # Create a FastaFile object of appropriate subclass.
    if options.format == "NCBI":
        fasta_in = NCBIFastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
    elif options.format == "UNIPROT":
        fasta_in = UniprotFastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
    elif options.format == "TIGR":
        fasta_in = TIGRFastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
    elif options.format == "MOUSE":
        fasta_in = MouseFastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
    elif options.format == "REFSEEK_CUSTOM":
        fasta_in = RefseekRNAFasta(fasta_file, experiment_id, options.source_db, options.protein_type)
    elif options.format == "PHYTOZOME":
        fasta_in = PhytozomeFastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
    else:
        fasta_in = FastaFile(fasta_file, experiment_id, options.source_db, options.protein_type)
        
    # Parse and import protein sequence records via FastaFile object.
    fasta_in.parse()
    fasta_in.db_import()

    # Clean up and exit
    print "Parsing fasta file {0} and importing to DB {1} complete.".format(fasta_file, engine.url)


if __name__ == "__main__":
    main()

