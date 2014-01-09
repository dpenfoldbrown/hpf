#!/usr/bin/python

# Script to translate the non-standard AA codes in already exported foldable sequences.
# Takes a directory of <pred_code>.fasta files, checks sequences in the files for non-
# standard codes, translates codes if found to create a new sequence, replaces file with 
# file containing translated sequence, and updates foldable record (hpf.filesystemOutfile) 
# in the hpf DB.
#
# dpb, 9/05/2012

import os
import re
import argparse

non_standard = {
    'X': 'A',   # Unknown -> Alanine
    'Z': 'Q',   # Glutamine/Glutamic Acid  -> Glutamine
    'B': 'N',   # Asparigine/Aspartic Acid -> Asparagine
    'U': 'C',   # Selenocysteine -> Cysteine
}
nonstandard_pattern = r"[{0}]".format("".join(non_standard.keys()))


def parse_code_from_file(file):
    """Parses the 2letter3number Rosetta prediction code from a filename"""
    file_pattern = r"(?P<code>[a-z]{2}[0-9]{3}).fasta"
    code_found = re.match(file_pattern, file)
    if not code_found:
        raise IOError("File '{0}' is not in exported foldable format".format(file))
    return code_found.group('code')

def parse_foldable_file(file):
    """Parses the sequence key and sequence from an exported foldable file"""
    header_pattern = r">sequence\.id\.(?P<sequence_key>[0-9]+).*"
    with open(file) as handle:
        header = handle.next()
        header_found = re.match(header_pattern, header)
        if not header_found:
            raise IOError("Header line '{0}...' in file {1} does not match foldable header format".format(header[:20], file))
        seq_key = header_found.group('sequence_key')
        sequence = handle.next().rstrip()
    return seq_key, sequence


def translate_foldables(dir, update_db=True):
    """Translates foldable fasta files in the manner described above. If update_db is true,
    adds new sequences translated to hpf.sequence table and updates hpf.filesystemOutfile 
    (foldable records) to link to the new sequence key.
    """
    from hashlib import sha1
    from hpf.hddb.reexport_foldables import write_fasta
    from hpf.hddb.db import push_to_db, Session, Sequence, FilesystemOutfile
   
    print "Translating foldable fastas found in dir {0}".format(dir)
    
    files = os.listdir(dir)
    for file in files:
        try: 
            code = parse_code_from_file(file)
        except IOError as e:
            print "{0}. Ignoring file..".format(e)
        sequence_key, sequence = parse_foldable_file(file)

        # If the sequence contains a non-standard AA code, translate nonstandard to normal codes
        if re.search(nonstandard_pattern, sequence):
            print "{0} contains nonstandard AA codes. Translating".format(file)
            print "\tOriginal  : {0}".format(sequence)
            
            translated_seq_id = "None"
            for nsaa in non_standard.keys():
                sequence = sequence.replace(nsaa, non_standard[nsaa])
            
            print "\tTranslated: {0}".format(sequence)
            
            if update_db:
                # Add new sequence to DB (push_ will return None if seq_dbo is already in DB)
                print "Adding translated sequence to the DB"
                session = Session()
                seq_dbo = Sequence(sequence=sequence, sha1=sha1(sequence).hexdigest())
                seq_dbo = push_to_db(session, seq_dbo, exception_str="Pushing sequence for code {0} failed".format(code), raise_on_duplicate=False)
                if not seq_dbo:
                    seq_dbo = session.query(Sequence).filter_by(sha1=sha1(sequence).hexdigest()).first()
                
                # Get foldable record and change seq key to new translated seq's id
                print "Updating foldable record from old ({0}) to new ({1}) sequence key".format(sequence_key, seq_dbo.id)
                foldable_dbo = session.query(FilesystemOutfile).filter_by(prediction_code=code).first()
                if not foldable_dbo:
                    raise Exception("Foldable record not found in DB for code {0}".format(code))
                foldable_dbo.sequence_key = seq_dbo.id
                session.flush()
                
                translated_seq_id = seq_dbo.id
            
            print "Writing translated foldable to file {0}".format(file)
            with open(file, 'w') as handle:
                write_fasta(handle, translated_seq_id, len(sequence), sequence)
    print "Translating foldables complete"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool to replace exported foldables with non-standard AA codes with standard code sequence files")
    parser.add_argument("-d", "--dir", action="store", dest="directory", required=True,
            help="Directory to read foldable fasta files from")
    parser.add_argument("--no_update", action="store_false", dest="update_db", default=True,
            help="Flag to NOT update the database with translated sequences. For testing/debugging")
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        raise EnvironmentError("Directory {0} is not a valid directory".format(args.directory))
    os.chdir(args.directory)

    translate_foldables(dir=args.directory, update_db=args.update_db)

