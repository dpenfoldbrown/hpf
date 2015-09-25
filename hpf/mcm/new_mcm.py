#!/usr/bin/env python

## MCM
# The whole of MCM functionality, derived from mcm.py, remcm.py, process.py, and mcm.pl
# Auth: dpb 3/28/2012
# Cluster setup functionality separate, in mcm_cluster_driver.py
# Ignoring task chain structure things, instead moving all functionality
# into appropriate modules, and running as a sequence of functions


import os
import re
import shutil
from hpf.hddb.db import Session, push_to_db
from hpf.pdb.psipred import HPFPsipredWrap, percent_alpha_beta
from hpf.mcm import McmDB, McmNoTaskFactory
from hpf.mcm.mammoth import MammothRun
from hpf.mcm.denovo_results import DenovoResultFile, DenovoResult
from hpf.mcm.cluster import RobettaClusterer, RobettaCluster, RobettaConvergence, RobettaClusterCenter

# RE for Rosetta denovo prediction code
CODE_PATTERN = r"(?P<code>[a-z]{2}[0-9]{3})"

# SCOP Version used
SCOP_VERSION = '1.75'

# Number of MCM records to keep (in DB) from each run. Keeps best.
KEEP_MCM = 5

class MCM():
    """A class for wrapping full MCM functionality. Completes all steps necessary for running MCM
    PARAMETERS:
    decoy_file  - Rosetta denovo results file (silent format), named by code. eg: na977.result
    code        - (Optional) Rosetta prediction code. Will be parsed from filename if not given.
                  NOTE: code given should match code in decoy filename
    work_dir    - (Optional) a base dir to create working directory in (work_dir/code)
    trim_score  - First filter for decoys from denovo results. Keep the best 'trim_score' scored decoys. DEFAULT 15,000
    trim_rg     - Second filter for denovo results decoys. Of the best 'trim_score' decoys, keep the best 'trim_rg' decoys. DEFAULT 10,000
                  NOTE: rg is Radius of Gyration. Lower is better.
    extract_log - A file for the pdb file extractor to output to. DEFAULT 'extract.log'
    rosetta_pathsfile - The location of the "paths.txt" file required for Rosetta (for extract functionality)
    mammoth_listfile  - The mammoth listfile containing all PDB filenames to compare against and their location
                          NOTE: Format is "MAMMOTH List\n<dir with pdb files>\n<pdbfile>\n<pdbfile>\n...<pdbfile>\n"
    mammoth_datafile  - The mammoth database file containing information on all astral scop structures
                          NOTE: Columns: "structure_key   sequence_key    length  percent_alpha   percent_beta    sccs    astral_ac"
    ginzu_version     - The ginzu_version of this MCM run. (The ginzu version of the results data that goes into this run)
                          NOTE: Used mostly for Psipred/Secondary structure prediction gathering
    nr_db       - The location of the blast 'nr' database. Used when Psipred SS preds do not exist for decoy sequences
                  in the HPF database.
    dbstore     - True stores cluster information (centers, convergence, and structures) and Mcm results in the database, False does not
    cleanup     - True removes working files, False leaves them in 'work_dir'/'code'
    debug       - True prints results information and debug information to screen

    INSTANCE VARIABLES:
    foldable_record - (dbo) ORM object of the foldable record corresponding to prediction code
    sequence        - (str) amino acid sequence OF THE FOLDABLE RECORD (note: may be different from domain sequence)
    sequence_key    - (int) sequence key OF THE FOLDABLE RECORD (see above), retrieved by Rosetta code (eg oa123)
    parent_sequence-key - (int) sequence key of the parent protein, from foldable record fetched by rosetta code
    """

    def __init__(self, decoy_file, 
                       code=None, 
                       work_dir='/tmp/mcm', 
                       trim_score=15000, 
                       trim_rg=10000, 
                       extract_log='extract.log',
                       rosetta_pathsfile='paths.txt', 
                       mammoth_listfile='list.mammoth', 
                       mammoth_datafile='data.mammothDb',
                       ginzu_version=4,
                       ignore_ginzu_version=False,
                       nr_db="nr",
                       dbstore=True, 
                       cleanup=True,
                       debug=False):
        #DEBUG
        print "Initializing total-MCM object...",
        
        self.decoy_file = os.path.abspath(os.path.expanduser(decoy_file))
        self.code = self._check_code(code) if code else self._parse_code(decoy_file)
        
        self.base_dir = os.path.abspath(os.path.expanduser(work_dir))
        self.work_dir = os.path.join(self.base_dir, self.code)

        self.trim_score = trim_score
        self.trim_rg = trim_rg
        self.extract_log = extract_log
        self.rosetta_pathsfile = rosetta_pathsfile

        self.mammoth_listfile = mammoth_listfile
        self.mammoth_datafile = mammoth_datafile
        self.ginzu_version = ginzu_version
        self.ignore_ginzu_version = ignore_ginzu_version 
        self.nr_db = nr_db

        self.dbstore = dbstore
        self.cleanup = cleanup
        self.debug   = debug
        
        # Define a session instance var for HPF DB use
        self.session = None
        
        # Retrieve foldable record, sequence, and sequence_key from HPF db via prediction code
        self.foldable_record = self._get_foldable_record()
        self.sequence = self.foldable_record.sequence.sequence
        self.sequence_key = self.foldable_record.sequence_key
        self.parent_sequence_key = self.foldable_record.parent_sequence_key

        #DEBUG
        print "Initializing total-MCM object complete"


    def run(self, ):
        """Instance variables set:
        filtered_decoy_file - filename of file that filtered set of decoys from original silent file is written to
        cluster_cmd - filename of the command file to be passed to clusterer
        cluster_out - filename of clusterer results file
        cluster_log - filename of log where clusterer STDOUT and STDERR are written
        ss_pred     - the Psipred secondary structure prediction string for this sequence
        """

        # Setup working environment (create dirs, etc)
        print "Setting up MCM working environment in {0}".format(self.work_dir)
        self._setup_env()
    
        # Parse decoy (denovo result) file and filter to X best scores, and Y of those best RGs
        print "Parsing and filtering de Novo results: {0} best score down to {1} best RG".format(self.trim_score, self.trim_rg)
        denovo_results = DenovoResultFile(filename=self.decoy_file, prediction_code=self.code)
        best_results = denovo_results.get_top_count(self.trim_score)
        best_results = sorted(best_results, key=lambda r: r.radius_gyration)[:self.trim_rg]

        # Write filtered results to file (to be passed to clusterer - the filtered decoy file will be used from here on)
        self.filtered_decoy_file = "{0}.score{1}.rg{2}".format(os.path.join(self.work_dir, self.code), self.trim_score, self.trim_rg)
        print "Writing filtered de Novo results to file {0}".format(self.filtered_decoy_file)
        denovo_results.write_to_file(outfile=self.filtered_decoy_file, results_list=best_results)

        # Parse filtered decoy file into a DenovoResultFile object and remove old denovo results object.
        # (Because clusterer runs on the filtered decoy file, need a DRF obj that corresponds to the indices output by the clusterer)
        filtered_denovo_results = DenovoResultFile(filename=self.filtered_decoy_file, prediction_code=self.code)
        del(denovo_results)

        # Run Rosetta clusterer on filtered set of decoys (not passing sequence - will be parsed from decoyfile)
        print "Running the robetta clusterer"
        self.cluster_cmd = self.code + ".cluster_cmd"
        self.cluster_out = self.code + ".cluster_out"
        self.cluster_log = self.code + ".cluster_log"
        clusterer = RobettaClusterer(command_file=self.cluster_cmd, 
                                     decoy_file=self.filtered_decoy_file, 
                                     outfile=self.cluster_out, 
                                     log_file=self.cluster_log, 
                                     path=self.work_dir)
        
        # Parse clusterer results into RobettaCluster object - contains RobettaConvergence obj (convergence) 
        #   and list+dict of RobettaClusterCenter objs (centers + centers_dict)
        print "Parsing clusterer results from file {0}".format(clusterer.outfile)
        cluster_results = RobettaClusterer.parse(clusterer.outfile)

        # Check returned cluster results: if no centers, problem in parsing
        if not cluster_results.centers:
            raise Exception("Error: no cluster centers parsed from outfile {0}".format(clusterer.outfile))

        #DEBUG: print cluster results
        if self.debug:
            print "+++Convergence:"
            print "\t{0}".format(cluster_results.convergence)
            print "+++Cluster Centers:"
            for c in cluster_results.centers:
                print "\t{0}".format(c)

        # Extract atom records from cluster centers and decoy file, output pdb files per center: decoy_<index>.pdb
        #   return is a dict, index => pdb_filename. Sets centers' pdb_file instance var
        print "Extracting cluster centers from de Novo silent result file into pdb atom record files"
        center_files_dict = self._extract_pdbs(filtered_denovo_results, cluster_results, self.extract_log, self.rosetta_pathsfile)

        # Store cluster info and centers in DB (rosetta_convergence (linked to fsO), rosetta_cluster, structure)
        if self.dbstore:
            print "Storing cluster convergence record and cluster centers to DB"
            convergence_record = self._store_convergence(cluster_results.convergence, self.foldable_record.id, self.filtered_decoy_file)
            center_record_dict = self._store_cluster_centers(cluster_results.centers, self.sequence_key, convergence_record.id)

        # Create hpf.mcm.McmDB object to contain mammoth list and data files, used in Mammoth and MCM functionality
        print "Creating Mammoth/MCM reference database object"
        with open(self.mammoth_datafile) as handle:
            mcmdb = McmDB.load(handle)
        mcmdb.list_file = self.mammoth_listfile
        mcmdb.scop = SCOP_VERSION

        # Mammoth cluster centers (via mcm.MammothRun) and parse mammoth scores (many per cluster center)
        print "Running Mammoth on all cluster centers, and parsing mammoth scores"
        mammoth = MammothRun(decoys=center_files_dict.values(), 
                             experiment_listfile=mcmdb.list_file, 
                             prediction_listfile="prediction.mammoth", 
                             directory=self.work_dir, 
                             outfile="mammoth.results", 
                             autorun=True,
                             debug=self.debug)
        mammoth_scores = mammoth.parse_scores()
        
        # Get foldable sequence's SS pred (get parent protein's SS, snip foldable's SS out of it) 
        #print "Getting Psipred SS predictions"
        fold_start, fold_stop = self._get_foldable_range(self.code, self.parent_sequence_key, self.sequence)
        self.ss_pred = self._get_ss_region(self.parent_sequence_key, self.ginzu_version, fold_start, fold_stop, ignore_ginzu_version=self.ignore_ginzu_version)
        if self.debug:
            print "Foldable range: {0} - {1}".format(fold_start, fold_stop)
            print "Foldable Psipred SS: {0}".format(self.ss_pred)
        
        # Calculate required MCM values (% alpha/beta, seq length) using the snipped-out SS prediction.
        sequence_length, percent_alpha, percent_beta = percent_alpha_beta(self.sequence, self.ss_pred)
        if self.debug: print "Percent alpha: {0}, Percent beta: {1}".format(percent_alpha, percent_beta)

        # MCM all mammoth scores (via McmFactory fctly). mcm_scores is list of McmData DBOs
        print "Creating McmFactory and running MCM on all mammoth scores"
        mcm_factory = McmNoTaskFactory(mcmdb, percent_alpha, percent_beta, cluster_results.convergence.radius1, debug=False)
        mcm_scores = map(mcm_factory.create, mammoth_scores)

        #DEBUG: print top MCM scores (number defined by KEEP_MCM)
        if self.debug:
            print "+++Top Kept MCM Scores (x2 to see top spread):"
            for m in sorted(mcm_scores, reverse=True)[:KEEP_MCM*2]:
                print "\t{0}".format(m.str_long())

        # Store mcm scores in DB (mcm/mcmdata table)
        if self.dbstore:
            print "Storing MCM records in HPF DB"
            # Sort in reverse: highest (aka best scores) first
            mcm_scores.sort(reverse=True)
            self._store_mcm_data(mcm_scores[:KEEP_MCM], self.sequence_key, self.foldable_record.id, convergence_record.id, center_record_dict)

        # Clean up working files
        if self.cleanup: 
            print "Removing working directory '{0}' and all contents".format(self.work_dir)
            self._cleanup()

        print "MCM run on {0}, sequence {1} complete".format(self.code, self.sequence_key)
        print "decoy file: {0}, filtered decoy file: {1}".format(self.decoy_file, self.filtered_decoy_file)
        print "working files in '{0}'".format(self.work_dir)


    def _parse_code(self, filename):
        """Parse Rosetta prediction code from given filename. If a code is not found, raise exception"""
        found = re.search(CODE_PATTERN, filename)
        if not found: 
            raise Exception("Could not parse prediction code from filename '{0}'".format(filename))
        return found.group('code')        

    def _check_code(self, code):
        """Checks code against static CODE PATTERN. If correct, returns code. Otherwise, exception"""
        if not re.match(CODE_PATTERN, code):
            raise Exception("Given code '{0}' does not match Rosetta prediction code form (eg: aa111)".format(code))
        return code

    def _setup_env(self, ):
        """Creates and checks directories and files. Changes working directory"""
        if not os.path.isdir(self.base_dir):
            os.mkdir(self.base_dir)
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        os.chdir(self.work_dir)

        # Check decoy file's existence
        if not os.path.isfile(self.decoy_file):
            raise IOError("Given decoy file '{0}' does not exist (is not a file)".format(self.decoy_file))
        # Check for Mammoth files (always required)
        if not os.path.isfile(self.mammoth_listfile):
            raise IOError("Given mammoth list file '{0}' does not exist or is not accessible".format(self.mammoth_listfile))
        if not os.path.isfile(self.mammoth_datafile):
            raise IOError("Given mammoth data file '{0}' does not exist or is not accessible".format(self.mammoth_datafile))
        # Check for rosetta_pathsfile, required for extracting pdb atom records from silent files
        if not os.path.isfile(self.rosetta_pathsfile):
            raise IOError("Given rosetta paths file '{0}' does not exist or is not accessible".format(self.rosetta_pathsfile))
	#kdrew: removing nr db check
        # Do a warning check for files required if Psipred will be run
        #if not os.path.isfile(self.nr_db):
        #    raise Warning("Given nr_db file '{0}' does not exist or is not accessible".format(self.nr_db))

    def _check_session(self, ):
        """Opens session if self.session is None. Closes the session (keeps it from expiring, auto-opened when used)"""
        if not self.session:
            self.session = Session()
        self.session.close()

    def _extract_pdbs(self, denovo_results, cluster_results, log=None, pathfile=None):
        """For each cluster center in cluster results, make an individual silent file and then extracts the PDB
        record from that silent file (outputs PDB record to file)
          denovo_results  - DenovoResultFile object (holds silent records)
          cluster_results - RobettaCluster object (holds cluster centers)
        Returns a dict of form {index => pdb_filename} for all cluster centers
        NOTE: It looks odd to create an individ. silent file for each decoy, but is the only way to guarantee
        extractor is pulling the right decoy for given S_id (which are not unique - why, I'll never know)
        """
        from hpf.mcm.extract import extract
        pdb_files = dict()
        
        for center in cluster_results:
            # Create filenames for holding silent file and target to move output PDB file to
            silent_file = "decoy_{0}.silent".format(center.index)
            pdb_file = "decoy_{0}.pdb".format(center.index)
            
            # Get and write silent record
            silent_record = denovo_results[center.index]
            denovo_results.write_to_file(outfile=silent_file, results_list=[silent_record])
            
            # Run extract and move extracted file to named pdb file
            extract_file = extract(silent_file, center.rosetta_id, log_file=log, paths=pathfile, debug=self.debug)
            shutil.move(extract_file, pdb_file)
            if not os.path.isfile(pdb_file):
                raise OSError("Extract functionality in MCM failed to create pdb file '{0}'".format(pdb_file))
           
            # Set center's pdb_file attribute to created file and add to dictionary
            center.pdb_file = pdb_file
            pdb_files[center.index] = pdb_file
            
            if self.cleanup: 
                os.remove(silent_file)
        return pdb_files
    
    def _get_foldable_record(self, ):
        """Gets an hpf.filesystemOutfile (ORM: FilesystemOutfile) record from the DB based on prediction code
        Should be considered a record of the foldable domain sequence (filesystemOutfile name makes no sense)
        """
        from hpf.hddb.db import FilesystemOutfile
        self._check_session()
        foldable = self.session.query(FilesystemOutfile).filter_by(prediction_code=self.code).first()
        if not foldable:
            raise Exception("Failed to find foldable record with code '{0}' in DB".format(self.code))
        return foldable
    
    def _get_foldable_range(self, prediction_code, parent_seq_key, foldable_sequence):
        """Returns the start and stop RANGE numbers (beginning at 1, not 0) of the foldable record's sequence
        (the position at which the foldable sequence starts and stops in the sequence of its parent protein). 
        Returns tuple of long ints if found
        """
        from hpf.hddb.db import Domain
        self._check_session()
        domain = self.session.query(Domain).filter_by(ibm_prediction_code=prediction_code, parent_sequence_key=parent_seq_key).first()
        if not domain:
            raise Exception("No domain found in DB with code {0} and parent seq key {1}".format(prediction_code, parent_seq_key))
        if not domain.region:
            raise Exception("Domain ID:{0} has no region. Database error".format(domain.id))
        
        f = re.search(foldable_sequence, domain.sequence.sequence)
        if not f:
            raise Exception("Foldable sequence '{0}' not found within Domain ID: {1} sequence".format(foldable_sequence, domain.id))
        fold_start_index = f.start()
        fold_stop_index = f.end() - 1
        return domain.region.start + fold_start_index, domain.region.start + fold_stop_index 
        
    def _get_ss_region(self, parent_sequence_key, ginzu_version, start, stop, ignore_ginzu_version=False):
        """Returns a region of the parent sequence's Psipred SS string. First checks
        the HPF DB for a pre-existing parent SS prediction. If found, return the appropriate 
        region. Otherwise, run Psipred on the parent sequence and snip the range out of
        the resulting SS string.
        NOTE: start and stop are "range" numbers. IE, 1-indexed (first residue at 1, second at 2...
        EG  : start 1, stop 10 will give substring starting at index 0, last char at index 9
        """
        from hpf.hddb.db import Psipred
        self._check_session()
        if ignore_ginzu_version:
            parent_psipred = self.session.query(Psipred).filter_by(sequence_key=parent_sequence_key).first()
        else:
            parent_psipred = self.session.query(Psipred).filter_by(sequence_key=parent_sequence_key, ginzu_version=ginzu_version).first()

        if parent_psipred:
            parent_ss = parent_psipred.prediction
        else:
            parent_ss = HPFPsipredWrap(sequence_key=parent_sequence_key,
                                          nr_db=self.nr_db,
                                          ginzu_version=self.ginzu_version,
                                          autorun=True,
                                          dbstore=self.dbstore
                                          ).get_prediction_string()
        #DEBUG
        if self.debug:
            print "Parent protein Psipred SS: {0}".format(parent_ss)

        return parent_ss[start-1:stop]


    def _store_convergence(self, convergence, foldable_key, decoy_file):
        """Stores convergence values in hpf.rosetta_convergence (ORM: RosettaConvergence)
        Links to hpf.filesystemOutfile via foldable_key (id)
        Parameters:
          convergence   - hpf.mcm.cluster.RobettaConvergence object containing cluster convergence info
          foldable_key  - the ID of the hpf.filesystemOutfile (ORM: FilesytemOutfile) entry for this code
          decoy_file    - the filename of decoy file given to the clusterer to cluster
        Returns the successfuly added RosettaConvergence ORM object
        """
        from hpf.hddb.db import RosettaConvergence
        self._check_session()
        cv = RosettaConvergence(outfile_key=foldable_key, target=decoy_file,
                radius1=convergence.radius1, size1=convergence.size1,
                radius2=convergence.radius2, size2=convergence.size2,
                total_decoys=convergence.total_decoys)
        cv = push_to_db(self.session, cv, exception_str="Failed to add RosettaConvergence (outfile_key: {0}) to DB".format(cv.outfile_key))
        return cv

    def _store_cluster_centers(self, cluster_centers, sequence_key, convergence_key):
        """Stores cluster centers in hpf.rosetta_cluster (ORM: RosettaCluster) and centers' structures in hpf.structure (ORM: Structure)
        Parameters:
          cluster_centers - list of hpf.mcm.cluster.RobettaClusterCenter objs
          sequence_key    - the sequence key of the domain sequence decoys were created from
          convergence_key - ID of the corresponding RosettaConvergence ORM object to link ClusterCenters to
        NOTE: RobettaCluster objects MUST have their pdb_file parameters set
        Returns dict of added RosettaCluster DBOs (linked to Structure objs by key), {index => RosettaCluster obj}
        """
        from hpf.hddb.db import Structure, RosettaCluster
        self._check_session()
        centers_dict = dict()
        for center in cluster_centers:
            atom_record = center.get_atom_record()
            if atom_record == None or atom_record == "":
                raise Exception("Failed to get atom record for center {0}".format(center))
            
            # Create and push Structure ORM object
            struct = Structure(sequence_key=sequence_key, structure_type="decoy", comment=center.pdb_file, text=atom_record)
            struct = push_to_db(self.session, struct, exception_str="Failed to add Structure for center {0} to DB".format(center))
            
            # Create and push RosettaCluster
            cc = RosettaCluster(index=center.index, size=center.size, rank=center.rank, convergence_key=convergence_key, structure_key=struct.id)
            push_to_db(self.session, cc, exception_str="Failed to add RosettaCluster for center {0} to DB".format(center))

            # Add to centers dict
            centers_dict[center.index] = cc
        return centers_dict

    def _store_mcm_data(self, mcm_scores, sequence_key, outfile_key, convergence_key, centers_dict):
        """Creates and populates Mammoth and corresponding McmData ORM objs with linking IDs (sequence, outfile,
        convergence, and structure), then pushes to the DB. Generally store top 5 MCM scores
        Parameters:
          mcm_scores    - list of hpf.hddb.db.McmData objects to store in DB
          sequence_key  - ID of sequence from which MCM scores come
          outfile_key   - ID of the foldable record (FilesystemOutfile) corresponding to the seq and MCM scores
          convergence_ke- ID of the convergence info for this sequence's cluster info (RosettaConvergence)
          centers_dict  - dict of form {CC index => OBJ} where object is anything w/ instance variable 'structure_key'
                          corresponding to the CC index's structure (EG: RosettaCluster DBO)
        This function will also fetch the structure key of the MCM score's mammothed cluster center and the
        structure key of the MCM score's mammoth astral structure (must do per mcm score)
        """
        from hpf.hddb.db import MammothFactory
        factory = MammothFactory()
        
        self._check_session()

        for score in mcm_scores:
            cc_index = int(score.mammoth.prediction.split(".")[0][6:])
            structure_key = centers_dict[cc_index].structure_key
            astral_structure_key = int(score.mammoth.experiment.split(".")[0])
            
            mammoth_dbo = factory.create(score.mammoth)
            mammoth_dbo.p_structure_key = structure_key
            mammoth_dbo.e_structure_key = astral_structure_key
            push_to_db(self.session, mammoth_dbo, 
                       exception_str="Failed to add {0} to DB for sequence {1}, index {2}".format(mammoth_dbo, sequence_key, cc_index))
            
            score.sequence_key = sequence_key
            score.outfile_key  = outfile_key
            score.convergence_key = convergence_key
            score.structure_key = structure_key
            score.astral_structure_key = astral_structure_key
            push_to_db(self.session, score,
                       exception_str="Failed to add {0} to DB for sequence {1}, index {2}".format(score, sequence_key, cc_index))
    
    def _cleanup(self, ):
        """Force removal of created working directory"""
        from subprocess import check_call
        if re.search(CODE_PATTERN+r"/?$", self.work_dir):
            ret = check_call(["rm", "-r", "-f", self.work_dir])
        else:
            raise Exception("Working directory '{0}' not valid for removal (must be a code directory, eg oa123".format(self.work_dir))
        return ret


