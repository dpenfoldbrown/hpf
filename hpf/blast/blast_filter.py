##
## blast_filter - A file containing classes for creating and running blast analyses and filters on
##   blast results. Contains methods for writing filtered records to files in Fasta format. Contains
##   subclasses to create more specific filters and filter results.
##
## Commentary - Duncan Penfold-Brown, 11/17/2010
##

## Imports

import re
from sqlalchemy import *
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio import SeqIO


## Classes

# Filtered:
# A class representing a query sequence and blast hit (match) of that sequence.
# This superclass stores only the ID of the query sequence (query_id) and the ID of a 
# single matching sequence (hit_id).
class Filtered():
    # __init__:
    # query_id - (str) The ID of the query sequence
    # hit_id   - (str) The ID of sequence that matches the query sequence
    def __init__(self,q_id, h_id):
        self.query_id = q_id
        self.hit_id = h_id

    # print_filt:
    # Prints a filtered object to stdout in the form QUERY_ID HIT_ID
    # Added to increase extensibility - accessing Filtered instance variables directly in Blast2SeqKey methods
    # increases coupling (changes here would require changes there).
    def print_filt(self):
        print self.query_id, self.hit_id

# BlastFilter:
# A class representing blast run and filter functionality. Contains instance variables for defining 
# a blast run and also parameters for filtering blast results. Contains functionality for outputting
# Fasta-formatted filter results to file.
class BlastFilter():
    
    def __init__(self, blast_exe=None, blast_db=None, blast_query=None, eval_cutoff = 1e-8, length_cutoff = .85, identity_cutoff = 0, blast_prog = "blastp", blast_processors=1, multi_hits=False):
        self.blast_exe = blast_exe
        self.blast_db = blast_db
        self.blast_query = blast_query
        self.blast_prog = blast_prog
        self.eval_cutoff = eval_cutoff
        self.length_cutoff = length_cutoff
        self.identity_cutoff = identity_cutoff
        self.blast_processors = blast_processors
        self.return_multiple_hits = multi_hits

    def runBlast(self, result_handle=None):
    # If a filehandle is given as input, simply reads and parses blast results from the input file into blast_records.
    # If a filehandle is not given as input, runs a new blast (with local arguments: blast_exe, blast_prog, etc.) on
    # Output: an iterator over a sequence of Record objects
        # If no filehandle given, or filehandle given is None, run new blast.
        if result_handle == None:
            result_handle, error_handle = NCBIStandalone.blastall(self.blast_exe, self.blast_prog, self.blast_db, self.blast_query, nprocessors=self.blast_processors)
        
        # Parse and return blast records from given filehandle or new blast run.
        blast_records = NCBIXML.parse(result_handle)
        return blast_records

    
    def filterBlast(self, blast_records):
    # Output: A dictionary of Filtered objects. Key: query_id of Filtered obj., Value: Filtered object or
    # list of Filtered Objects (all with the same query_id).
        
        filter_dict= {}
        num_given    = 0
        num_filtered = 0

        # For each record in blast_records, try to filter the record. If record passes filter, add to filter_dict.
        for blast_record in blast_records:
            num_given += 1
            filtered_record = self.filterBlastOne(blast_record)
            if filtered_record != None:
                # If multiple matches are desired, add the list of Filtered objects returned by filterBlastOne to the dict.
                # The key for a list of Filtered objects is the query_id of the first list object (query IDs the same for multiple matches).
                num_filtered += 1
                if self.return_multiple_hits:
                    #DEBUG
                    #print "blast_filter:: Query {0} returned {1} filtered records".format(blast_record.query, len(filtered_record))
                    filter_dict[str(filtered_record[0].query_id)] = filtered_record
                else:
                    filter_dict[str(filtered_record.query_id)] = filtered_record
            #else:
                #DEBUG
                #print "blast_filter:: Query {0} returned NO filtered records".format(blast_record.query)

        #DEBUG
        print "blast_Filter:: Number records given to filter: {0}".format(num_given)
        print "blast_filter:: Number records after filter: {0}".format(num_filtered)

        # Return the dictionary of filter results.
        return filter_dict
    
    def parse_blast_hit(self, record_query, align_def):
    # Parses input query ID and hit ID values. Parses values, then creates a Filtered object
        query_id = record_query.split('|')[0].strip()
        hit_id = align_def.split()[0].strip()

        return Filtered(query_id, hit_id)

    # filterBlastOne:
    # Applies a filter specified by instance variable parameters to the passed in blast record.
    # If the record meets filter requirements, a Filtered object of the record (or a list of Filtered 
    # objects, if return_multiple_hits is true) is returned. Filtered objects are created by calling
    # the parse_blast_hit local method.
    # Input: blast_record - A blast record of type Record (via Bio.Blast).
    # Ouput: A Filtered object or list of Filtered objects that meet filter requirements.
    # Changed: Loop structure changed and variables removed to increase loop efficiency (optimized for
    # returning single hits).
    def filterBlastOne(self, blast_record):
        # Create an empty list to store Filtered objects representing matches that meet filter requirements.
        filter_rec_list = list()
        
        # For each alignment in the blast record, if the alignment meets filter requirements, store it.
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.eval_cutoff and hsp.align_length > alignment.length*self.length_cutoff \
                    and hsp.identities > self.identity_cutoff:
                    
                    # If only a single hit is desired per query, return first match immediately.
                    if not self.return_multiple_hits:
                        return self.parse_blast_hit(blast_record.query, alignment.hit_def)
                    # Otherwise, append hit to the filter_rec_list and continue
                    else:
                        filter_rec_list.append(self.parse_blast_hit(blast_record.query, alignment.hit_def))

        # Return list of Filtered objects representing blast hits, or None (if nothing meets filter reqs)
        if len(filter_rec_list) > 0:
            return filter_rec_list
        else:
            return None

    # print_filtered_records:
    # Prints a dictionary of Filtered objects.
    # Input: dictionary of the form {(str)query_id: Filtered object} OR {(str)query_id: List of Filtered objects}.
    # Output: print to stdout.
    def print_filtered_records(self, filt_dict):
        for filt in filt_dict.values():            
            # If the dictionary value is a list of Filtered objects, print all matches in list.
            if isinstance(filt, list):
                for match in filt:
                    match.print_filt()
            else:
                filt.print_filt()

    # writeFasta:
    # Finds the sequence records in the local blast_query file that appear in the input filtered record
    # dictionary, and writes these Fasta-formatted sequence records into a given (or generated) outfile.
    # Input: filter_dict - A dictionary of Filtered objects.
    #        outfile (Optional) - a filehandle to store filtered sequences in Fasta format.
    def writeFasta(self, filter_dict, out_file=None):
        #print self.blast_query
        #print filter_dict
        
        # Open the blast_query file (given in object construction).
        handle_in = open(self.blast_query)
        
        # Open the outfile for writing. If one isn't provided, create one ('blast_query'_filtered.fasta).
        if None == out_file:
            handle_out = open(self.blast_query.rpartition('.')[0]+"_filtered.fasta", "w")
        else:
            handle_out = open(out_file, "w")
            
        # Create an empty list to store sequence records that appear in the input filter_dict.
        out_list = []
        
        # For every sequence record in the blast_query file, if the sequence record is in the filter_dict,
        # add the sequence record to out_list
        for seq_record in SeqIO.parse(handle_in, "fasta") :
            #print seq_record.id
            #print repr(seq_record.seq)
            #print len(seq_record)
            if seq_record.id.split('|')[0] in filter_dict:
                print seq_record.id
                out_list.append(seq_record)
        
        # Write all sequence records in outlist (records that are in both the blast_query file and filtered_dict) to outfile.
        SeqIO.write(out_list, handle_out, "fasta")
        
        # Close working files.
        handle_in.close()
        handle_out.close()


class DetailedFiltered(Filtered):
# DetailedFiltered:
# Extends Filtered class to support evalue and bit-score instance variables.

    def __init__(self, q_id, h_id, h_eval, h_bitscore):
        Filtered.__init__(self, q_id, h_id)
        self.hit_evalue = h_eval
        self.hit_bitscore = h_bitscore
    
    # Prints all the information stored in a DetailedFiltered object.
    def print_full(self):
        print self.query_id, self.hit_id, self.hit_evalue, self.hit_bitscore

class DetailedBlastFilter(BlastFilter):
# Extends BlastFilter class to support better formatting and capturing of blast hit information
# (via regular expression use) in parse_blast_hit. Also extends functionality for retrieving evalue
# and bit-score fields for each query match that meets filter requirements (in filterBlastOne).
# Also works with DetailedFiltered objects, an extension of Filtered objects.

    def _threshold(self, eval, align_percent, identity):
    # A threshold function. Returns true if e-value, alignment percent, and identity cutoffs are met
    # Meant to match the phrase "take matches with id_cut identity over align_percent over the sequence"
        return float(eval) <= float(self.eval_cutoff) and \
               float(align_percent) >= float(self.length_cutoff) and \
               float(identity) >= float(self.identity_cutoff)
    
    def filterBlastOne(self, blast_record):
        
        # For each alignment in the blast record, if the alignment meets filter requirements, store it.
        filter_rec_list = list()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                
                align_length = max(len(hsp.query), len(hsp.match))
                max_seq_len  = max(blast_record.query_letters, alignment.length)
                align_percent = float(align_length) / float(max_seq_len)
                identity_percent = float(hsp.identities) / float(blast_record.query_letters)
                
                if self._threshold(hsp.expect, align_percent, identity_percent):
                    if not self.return_multiple_hits: 
                        return self.parse_blast_hit(blast_record.query, alignment.hit_def, hsp.expect, hsp.bits)
                    else:
                        filter_rec_list.append(self.parse_blast_hit(blast_record.query, alignment.hit_def, hsp.expect, hsp.bits))

        # Return list of Filtered objects representing blast hits, or None (if nothing meets filter reqs)
        if len(filter_rec_list) > 0:
            return filter_rec_list
        else:
            return None
    
    # Supports different parsing of align_def value to retrieve hit_id. Also creates a DetailedFiltered
    # object to support more alignment information (for eventual DB population).
    def parse_blast_hit(self, record_query, align_def, e_val, bit_score):
        # Query ID string: 'YER177W'. Can be captured as is.
        # Hit ID string: 'bddb|273260| Saccharomyces cerevisiae [BDDB-SEQID 273260]'. Parse and capture. 
        
        query_id = record_query.strip()
        hit_re = "bddb\|(?P<hit_id>.*?)\|.*$"
        match = re.search(hit_re, align_def)
        hit_id = match.group('hit_id')
        
        # Return a DetailedFiltered object (NOTE: e-value and bit-score are floats).
        return DetailedFiltered(query_id, hit_id, e_val, bit_score)

    # write_filtered_DB:
    # Writes the contents of a dictionary of DetailedFiltered objects into a mySQL database specified by
    # the input db_string (otherwise, to localhost/test).
    # Specifically creates a table (input, defaults to alignments) containing rows of the form
    #   (query_id,  hit_id,  e-value,  bit-score,) 
    # where (query_id, hit_id) is a composite primary key. Rows are populated from DetailedFiltered objects. 
    # NOTE: WILL OVERWRITE TABLES AN EXISTING TABLE OF THE SAME NAME.
    # Input: filt-dict - Dictionary of the form {(str)query_id: Filtered object OR List of Filtered objects}.
    #   db_string - (str) String of the form 'username:password@host:port/database' specifying the DB to use.
    #   table - (str) The name of the table to populate with filter alignments.
    #   db_echo - (bool) A switch for printing DB control statements to screen. False: no prints.
    # Output: To given or default database.
    def write_filtered_DB(self, filt_dict, db_string, table, db_echo = False):
        print "Populating MySQL database {0}, table {1} with sequence alignment information.".format(db_string, table)
        
        # Create MySQL database engine (responsible for communication with db service).
        db_engine = create_engine('mysql+mysqldb://'+db_string, echo = db_echo)
        
        # Create metadata catalog for holding information on the database, and bind the engine to the metadata (for ease of use).
        db_meta = MetaData()
        db_meta.bind = db_engine
        
        # Create table for storing query and filtered hit information. Will contain a field for every instance variable of DetailedFiltered.
        alignments = Table(table, db_meta, 
            Column('query_id', String(100), primary_key = True),
            Column('sequence_key', String(100), primary_key = True),
            Column('e_value', Float),
            Column('bit_score', Float),
            )
        # Note that this statement will overwrite an existing table of the same name in the given database.
        alignments.drop(db_engine, checkfirst = True)
        alignments.create(db_engine)
        
        # Insert all data from all objects in given dictionary into database
        for filt in filt_dict.values():            
            # If the dictionary value is a list of DetailedFiltered objects, process all matches in list.
            if isinstance(filt, list):
                for match in filt:
                    alignments.insert().execute(query_id = match.query_id, sequence_key = match.hit_id, \
                        e_value = match.hit_evalue, bit_score = match.hit_bitscore)
            else:
                alignments.insert().execute(query_id = filt.query_id, sequence_key = filt.hit_id, \
                    e_value = filt.hit_evalue, bit_score = filt.hit_bitscore)


# AstralFiltered:
# Extends Filtered class to support another instance variable, superfamily, to be given on instantiation.
class AstralFiltered(Filtered):
    def __init__(self,q_id,h_id,h_sf):
        Filtered.__init__(self,q_id,h_id)
        self.superfamily = h_sf

# AstralBlastFilter:
# Extends BlastFilter class to support correct parsing of blast records in order to retrieve superfamily value.
# Causes output of blastFilter to be a dictionary of AstralFiltered objects or lists of AstralFiltered objects.
class AstralBlastFilter(BlastFilter):
    def parse_blast_hit(self, record_query, align_def):
        query_id = record_query.split('|')[0].strip()
        hit_id = align_def.split()[0].strip()
        sccs_id = align_def.split()[1].strip()
        superfamily = sccs_id.rpartition('.')[0]

        return AstralFiltered(query_id, hit_id, superfamily)

# GOBlastFilter:
# Extends BlastFilter class to print raw query_id and hit_id fields before parsing and creating Filtered objects.
class GOBlastFilter(BlastFilter):
    def parse_blast_hit(self, record_query, align_def):
        print "query_id_raw: ", record_query
        query_id = record_query.split('|')[1].strip()
        print "hit_id_raw: ", align_def
        hit_id = align_def.split('|')[0].strip()
        print query_id, hit_id

        return Filtered(int(query_id), int(hit_id))

class UniprotBlastFilter(DetailedBlastFilter):
# Parses result fields from standard uniprot headers

    def parse_blast_hit(self, record_query, hit_def, e_val, bit_score):
        # Query ID string: 'KYOTOGRAIL2005.734.1.1|MKG.734.1.s0_at'
        # Hit ID string  : 'tr|F6ZD22|F6ZD22_CIOIN Uncharacterized protein (Fragment) OS=Ciona intestinalis GN=Cin.39466 PE=4 SV=1' 
        query_id = record_query
        
        hit_pattern = r"[a-zA-Z]+\|(?P<hit_id>[a-zA-Z0-9-_]+)\|.*"
        hit_found = re.match(hit_pattern, hit_def)
        if not hit_found:
            raise Exception("Match header '{0}' of unknown format".format(hit_def))
        hit_id = hit_found.group('hit_id')
        
        # Return a LocalHumanFiltered object
        return DetailedFiltered(query_id, hit_id, e_val, bit_score)



class LocalHumanFiltered(DetailedFiltered):

    def __init__(self, query_id, hit_id, hit_protein_id, hit_experiment_id, hit_eval, hit_bitscore):
    # With superclasses, field will be:
    # self. query_id, hit_id, hit_evalue, hit_bitscore, hit_protein_id, hit_experiment_id
        DetailedFiltered.__init__(self, query_id, hit_id, hit_eval, hit_bitscore)
        self.hit_protein_id = hit_protein_id
        self.hit_experiment_id = hit_experiment_id
    
class LocalHumanBlastFilter(DetailedBlastFilter):
# Where query is Kevin's local export format, results are dpb's local export format

    def parse_blast_hit(self, record_query, hit_def, e_val, bit_score):
        # Query ID string: 'bddb|4225364| [BDDB-SEQID 4225364]' or 'bddb|5167483| [BDDB-SEQID 5167483 - Human RNA extra]'
        # Hit ID string  : '12006 | 316230 | 804 | KnownProteins EMBL Release 22.34d.1 (H. sapiens)' 
        
        query_pattern = r"bddb\|(?P<query_id>[0-9]+)\|.+"
        hit_pattern = r"(?P<hit_id>[0-9]+) \| (?P<protein_id>[0-9]+) \| (?P<experiment_id>[0-9]+) .*"

        query_found = re.match(query_pattern, record_query)
        if not query_found:
            raise Exception("Query header '{0}' of unknown format".format(record_query))
        query_id = query_found.group('query_id')

        hit_found = re.match(hit_pattern, hit_def)
        if not hit_found:
            raise Exception("Match header '{0}' of unknown format".format(hit_def))
        hit_id = hit_found.group('hit_id')
        hit_protein_id = hit_found.group('protein_id')
        hit_experiment_id = hit_found.group('experiment_id')
        
        # Return a LocalHumanFiltered object
        return LocalHumanFiltered(query_id, hit_id, hit_protein_id, hit_experiment_id, e_val, bit_score)

class HumanRNABlastFilter(DetailedBlastFilter):
# For the landthaler Human RNA data from berlinsies.

    def parse_blast_hit(self, record_query, hit_def, e_val, bit_score):
        # Query ID string: 'NP_002810.1 polypyrimidine tract-binding protein 1 isoform a [Homo sapiens].'
        # Hit ID string  : '12006 | 316230 | 804 | KnownProteins EMBL Release 22.34d.1 (H. sapiens)' 
        
        query_pattern = r"(?P<query_id>[A-Z]{2}_[0-9]+\.[0-9]+)\s.*"
        hit_pattern = r"(?P<hit_id>[0-9]+) \| (?P<protein_id>[0-9]+) \| (?P<experiment_id>[0-9]+) .*"

        query_found = re.match(query_pattern, record_query)
        if not query_found:
            raise Exception("Query header '{0}' of unknown format".format(record_query))
        query_id = query_found.group('query_id')

        hit_found = re.match(hit_pattern, hit_def)
        if not hit_found:
            raise Exception("Match header '{0}' of unknown format".format(hit_def))
        hit_id = hit_found.group('hit_id')
        hit_protein_id = hit_found.group('protein_id')
        hit_experiment_id = hit_found.group('experiment_id')
        
        # Return a LocalHumanFiltered object
        return LocalHumanFiltered(query_id, hit_id, hit_protein_id, hit_experiment_id, e_val, bit_score)

class HPDFiltered(DetailedFiltered):

    def __init__(self, query_id, query_start, query_stop, hit_id, hit_protein_seqkey, hit_experiment_id, hit_dom_type, hit_eval, hit_bitscore):
    # With superclasses, field will be:
    # self. query_id, hit_id, hit_evalue, hit_bitscore
        DetailedFiltered.__init__(self, query_id, hit_id, hit_eval, hit_bitscore)
        self.query_start = query_start
        self.query_stop = query_stop
        self.hit_protein_seqkey = hit_protein_seqkey
        self.hit_experiment_id = hit_experiment_id
        self.hit_domain_type = hit_dom_type

class HPDBlastFilter(DetailedBlastFilter):
# For the HPD list of domains from Jamboree
    def parse_blast_hit(self, record_query, hit_def, e_val, bit_score):
        # Query header: 'sp|A0AUZ9|CB067_HUMAN Uncharacterized protein C2orf67 OS=Homo sapiens GN=C2orf67 PE=1 SV=2:_: Range:142-330'
        # Hit header: '8374558|15196|804|msa<'

        query_pattern = r"[a-zA-Z]+\|(?P<uniprot>[a-zA-Z0-9_-]+)\|.*Range:(?P<start>[0-9]+)-(?P<stop>[0-9]+)"
        hit_pattern = r"(?P<dom_seq_key>[0-9]+)\|(?P<prot_seq_key>[0-9]+)\|(?P<exp_key>[0-9]+)\|(?P<dom_type>[a-zA-Z_-]+)"

        query_found = re.match(query_pattern, record_query)
        if not query_found:
            raise Exception("Query header '{0}' of unknown format".format(record_query))
        query_id = query_found.group('uniprot')
        query_start = query_found.group('start')
        query_stop = query_found.group('stop')
        query_id = "{0}:{1}-{2}".format(query_id, query_start, query_stop)

        hit_found = re.match(hit_pattern, hit_def)
        if not hit_found:
            raise Exception("Match header '{0}' of unknown format".format(hit_def))
        hit_id = hit_found.group('dom_seq_key')
        hit_protein_seqkey = hit_found.group('prot_seq_key')
        hit_experiment_id = hit_found.group('exp_key')
        hit_dom_type = hit_found.group('dom_type')

        # Return a LocalHumanFiltered object
        return HPDFiltered(query_id, query_start, query_stop, hit_id, hit_protein_seqkey, hit_experiment_id, hit_dom_type, e_val, bit_score)



