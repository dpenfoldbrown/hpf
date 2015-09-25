## A module of classes for handling HHPred/search results
# dpb 5/04/2012

import re
import os

class HHPredResultFile(object):
    """ A class to hold an HHPred (hhsearch) results file. 
    For explanation of header and hit fields, see HHPred output format documentation
    For more examples of database strings, see http://toolkit.tuebingen.mpg.de/hhpred
    Variables:
        [all header fields stored in instance variables]
        hits    - a list of HHPredHit objs representing hits in results file
    """

    def __init__(self, results_file, database, autoparse=True, debug=True):
        """ Parameters:
        results_file    - the filename of the results file. Must be openable
        database        - string describing the database the query was run against (eg scop175, scop70_175, pdb70_5Apr12...)
        autoparse       - will automatically parse hits and populate self.hits list
        """
        self.results_file = os.path.abspath(os.path.expanduser(results_file))
        self.database = database
        self.debug = debug

        # Results file header fields (populated in _parse_header())
        self.query = None
        self.match_columns = None
        self.num_seqs = None
        self.neff = None
        self.hmms_searched = None
        self.date = None
        self.command = None
        self._parse_header()

        self.hits = list()
        if autoparse:
            self.hits = HHPredResultFile.parse(self.results_file, parent=self, debug=self.debug)

    def __repr__(self, ):
        return "<HHPredResultFile: {0} against {1}>".format(self.query, self.database)

    def _parse_header(self, ):
        """ Open results file and parse header information, store in instance vars. Stop when the
        Hits column header is reached
        """
        query_pattern = r"Query\s+(?P<query>.+)"
        numseqs_pattern = r"No_of_seqs\s+(?P<num_seqs>.+)"
        date_pattern = r"Date\s+(?P<date>.+)"
        command_pattern = r"Command\s+(?P<command>.+)"
        stop_pattern = r"\s+No\s+Hit\s+Prob\s+E-value\s+P-value\s+Score"

        with open(self.results_file) as handle:
            for line in handle:
                line = line.rstrip()
                if line.startswith('Query'):
                    self.query = re.match(query_pattern, line).group('query')
                elif line.startswith('Match_columns'):
                    self.match_columns = int(line.split()[1])
                elif line.startswith('No_of_seqs'):
                    self.num_seqs = re.match(numseqs_pattern, line).group('num_seqs')
                elif line.startswith('Neff'):
                    self.neff = line.split()[1]
                elif line.startswith('Searched_HMMs'):
                    self.hmms_searched = int(line.split()[1])
                elif line.startswith('Date'):
                    self.date = re.match(date_pattern, line).group('date')
                elif line.startswith('Command'):
                    self.command = re.match(command_pattern, line).group('command')
                elif line == "":
                    continue
                elif re.match(stop_pattern, line):
                    if self.debug: print "Parsing {0} header complete".format(self.results_file)
                    break
                else:
                    if self.debug: print "Line '{0}' not recognized. Skipping..".format(line)
    
    @staticmethod
    def parse(results_file, parent=None, debug=True):
        """Parses the given HHPred results file. Returns a list of HHPredHit objects
        parent  - HHPredResultFile object representing the results_file (source of HHPredHit objs)
                  Passed in as opposed to self referenced to support static parsing
        debug   - a debug flag, default True.
        """
        hits = list()
        hit_pattern = r"^\s*[0-9]+.*\([0-9]+\)$"
        
        with open(results_file) as handle:
            for line in handle:
                # Break when first alignment line encountered
                if line.startswith('>'):
                    if debug: print "Parsing {0} hits complete. Hits: {1}".format(results_file, len(hits))
                    break
                # Create and save hit obj when hit line matched
                elif re.match(hit_pattern, line):
                    hits.append( HHPredHit.from_string(line, parent) )
        return hits

    def sort_by_probability(self, ):
        """Sort hits by probabilty (default for HHPredHit obj)"""
        self.hits.sort()
        return self.hits
    
    def sort_by_eval(self, ):
        """Sort hits by evalue"""
        self.hits.sort(key=lambda h: h.evalue)
        return self.hits
    
    def sort_by_number(self, ):
        """Sort hits by their hit number"""
        self.hits.sort(key=lambda h: h.number)
        return self.hits


class HpfHHPredResultFile(HHPredResultFile):
    """A subclass of HHPredResultFile that requires an hhpred version key and parses sequence_key from file
    Creates a 'dbo' variable (as a property) to represent the corresponding hpf.hddb.db ORM object from the DB.
        'dbo' is None if HHPredResultFile is not represented in DB, otherwise is the corresponding ORM object.
    Adds the 'db_store' convenience method to push this object to the DB.
    Adds session instance variables to avoid redundant session creation
    """
    
    def __init__(self, version_key=None, *args, **kwargs):
        if not version_key:
            raise Exception("Hpf HHPredResultFile must be provided with a version key. See hhpred_version DB table")
        super(HpfHHPredResultFile,self).__init__(*args, **kwargs)
        self.version_key = version_key
        self.sequence_key = self._parse_seqkey(self.query)
        self.session = None
        self._dbo = None

    def _get_dbo(self, ):
        if not self._dbo:
            self._dbo = self._fetch_dbo()
        return self._dbo
    def _set_dbo(self, value):
        self._dbo = value
    dbo = property(_get_dbo, _set_dbo)
    
    def _parse_seqkey(self, str):
        """Attempts to parse an HPF sequence key from given string. This is set to parse 
        HHPred 'query' field of form 'hpf|<seqkey>|<proteinkey>|<experimentkey>|<experiment>'
        """
        return int( str.split('|')[1] )

    def _setup_db(self, ):
        if not self.session:
            from hpf.hddb.db import ScopedSession
            self.session = ScopedSession()
    
    def _fetch_dbo(self, ):
        """Queries the database for an already existing HHPRF (same sequence key and version).
        If query find matching HHPRF entry, returns it. If no matching entries in DB, returns None.
        """
        from hpf.hddb.db import HHPredResultFile as HHPRF
        if not self.session:
            self._setup_db()
        rf = self.session.query(HHPRF).filter_by(sequence_key=self.sequence_key, version_key=self.version_key).first()
        return rf

    def db_store(self, ):
        """Pushes the HHPredResultFile object to the HPF database (via hpd.hddb.db.HHPredRFFactory)
        Flow: get session; if HHPredRF object NOT already represented in DB: make ORM object; push to DB;
        set dbo property to pushed object; return dbo object
        Returns the added ORM object
        """
        if not self.session:
            self._setup_db()
        if self.dbo:
            print "HHPredResultFile already represented in DB: {0}. Returning...".format(self.dbo)
        else:
            from hpf.hddb.db import HHPredRFFactory, push_to_db
            hrf_dbo = HHPredRFFactory().create(self, debug=self.debug)
            push_to_db(self.session, hrf_dbo, exception_str="Failed to add {0} to the DB".format(hrf_dbo))
            self.dbo = hrf_dbo
        return self.dbo

    def db_store_hits(self, ):
        """A convenience function to push all self.hits Hit objects to the DB via HHPredHit.db_store()
        Raises warning if no hits are available.
        Raises warning if # of hits parsed (in self.hits) is different from # of hits in DB
        If HHpredResultFile object has not yet been pushed to DB, pushes it.
        Returns list of HHPredHit ORM objects.
        """
        if not self.hits:
            raise Warning("No hit objects to push to DB")
        if not self.dbo:
            self.db_store()
        if not self.session:
            self._setup_db()
        rf_key = self.dbo.id
        hit_dbos = list()
        for hit in self.hits:
            hit_dbos.append( hit.db_store(self.session, rf_key) )
        if self.debug:
            print "Hits parsed: {0}, Hits in DB: {1}".format(len(self.hits), len(hit_dbos))
        if len(self.hits) != len(hit_dbos):
            raise Warning("Number of hits parsed from file does not match number added to DB")
        return hit_dbos


class HHPredHit(object):
    """ A class to represent a single HHPred hit.
    For explanation of hit fields, see HHPred/HHSearch output format documentation.
    Subclass this in order to define custom 'hit' parsing method - can support any hit 
    string format, parsing any needed ID and description values from the hit string.
    'db_store' method is a convenience method for creating an ORM object and pushing
        it to the HPF DB
    """
    
    def __init__(self, parent_result, number, hit, probability, evalue, pvalue, score, ss, columns, query_range, template_range, match_states):
        """ Parameters are as found in HHPred results files (same order):
          parent_result   - HHPredResult obj. The result object from which this Hit was parsed.
          number  - int. Index number of hit. HHPred sorts by probability, so 1 will be highest probablity
          hit     - str. HHPred's hit field - a 30-char description of the hit (contains IDs - subclass to change parsing)
          probability     - float. 
          evalue         - float.
          pvalue         - float.
          score           - float. HHPred score. Does not correspond in ordering to probability
          ss              - float. HHPred's secondary structure score for the query
          columns         - int.
          query_range     - str.
          template_range  - str.
          num_match_states- int.
        Additional Variables:
          hit_id          - str. An ID for the hit parsed from the whole hit line
          hit_description - str. The remaining description from the hit field
          query_start     - int. The start value parsed from query_range
          query_stop      - int. ""
          template_start  - int. The start value parsed from template_range
          template_stop   - int.
        """
        self.parent_result = parent_result
        self.number = int(number)
        self.hit = hit
        self.probability = float(probability)
        self.evalue = float(evalue)
        self.pvalue = float(pvalue)
        self.score = float(score)
        self.ss = float(ss)
        self.columns = int(columns)
        self.query_range = query_range
        self.template_range = template_range
        self.match_states = int(match_states)

        self.hit_id, self.hit_description = self._parse_hit(self.hit)
        self.query_start, self.query_stop = self._parse_range(self.query_range)
        self.template_start, self.template_stop = self._parse_range(self.template_range)

    def __repr__(self, ):
        return "<HHPredHit: {0}, {1} (query {2}-{3}, template {4}-{5})\t{6}\t{7}>".format(self.number, self.hit_id, 
                    self.query_start, self.query_stop, self.template_start, self.template_stop, 
                    self.probability, self.evalue)

    def __cmp__(self, other):
        """Sort by probabilty, highest first"""
        return cmp(other.probability, self.probability)

    def _parse_hit(self, hit):
        """Parses and returns the hit_id and hit_description from a HHPred result line hit field
        NOTE: Subclass this class and override THIS METHOD to define custom hit_id/description schemes
        This class works for SCOP hit descriptions. ID is the astral structure ID. Description is
        the remaining hit field
        """
        # TODO: This could be more sophisticated and actually check for an ID...
        # SCOP hit: 'd1yksa2 c.37.1.14 (A:325-623) '
        # SCOP hit: 'd1q0ua_ c.37.1.19 (A:) Probabl'
        hit_id, hit_desc = hit.split(" ", 1)
        return hit_id, hit_desc

    def _parse_range(self, range):
        """Parses the int start and stop values from a given range, form <start>-<stop>"""
        range_pattern = r"(?P<start>[0-9]+)-(?P<stop>[0-9]+)"
        range_found = re.search(range_pattern, range)
        if not range_found:
            raise Exception("Range not found in string '{0}'".format(range))
        return int(range_found.group('start')), int(range_found.group('stop'))

    @staticmethod
    def from_string(hit_str, parent_result):
        """Creates and returns an HHPredHit object from a complete HHPred results line string
        parent_result can be set to None (if keeping track of parent HHPredResultFile unimportant)
        """
        # String: No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
        # String:  1 d1yksa2 c.37.1.14 (A:325-623)  100.0 1.1E-44       0  330.8 -10.3  259  738-1040   34-299 (299)
        hit_str = hit_str.rstrip()
        
        float_pat = r"[-+.0-9eE]+"
        hit_pattern  = r"\s*(?P<num>[0-9]+)\s+(?P<hit>.{30})\s+"
        hit_pattern += r"(?P<prob>{0})\s+(?P<eval>{0})\s+(?P<pval>{0})\s+(?P<score>{0})\s+".format(float_pat)
        hit_pattern += r"(?P<ss>{0})\s+(?P<cols>[0-9]+)\s+".format(float_pat)
        hit_pattern += r"(?P<qrange>[-0-9]+)\s+(?P<trange>[-0-9]+)\s*\((?P<matchstates>[0-9]+)\)"

        # Groups in hit RE: num, hit, prob, eval, pval, score, ss, cols, qrange, trange, matchstates
        found = re.match(hit_pattern, hit_str)
        if not found:
            raise Exception("Result string '{0}' does not match result pattern".format(hit_str))
        fields = list(found.groups())

        # Using list expansion on RE capture list to instantiate
        return HHPredHit(parent_result, *fields)

    def db_store(self, session, resultfile_key):
        """A convenience method to create HHPredHit ORM objects and push them to the HPF DB
        given a pre-existing session. Most easily called with a HHPredResultFile parent object
        """
        from sqlalchemy.exc import IntegrityError
        from hpf.hddb.db import HHPredHitFactory, HHPredHit as HHPH, push_to_db
        hit_dbo = HHPredHitFactory().create(self, resultfile_key=resultfile_key)
        try:
            push_to_db(session, hit_dbo, exception_str="Failed to add HHPredHit {0} to the DB".format(hit_dbo))
        except IntegrityError:
            print "HHpredHit {0} already exists in DB. Returning pre-existing ORM object".format(hit_dbo)
            hit_dbo = session.query(HHPH).filter_by(resultfile_key=hit_dbo.resultfile_key, number=hit_dbo.number).first()
        return hit_dbo
        
