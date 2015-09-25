#!/usr/bin/env python

# Module file containing interpro class to run interpro on the commandline

import os
import subprocess

IPRSCAN = "iprscan"

class InterproResultFile():

    def __init__(self, filename, out_format="raw", iprlookup=True, goterms=True, autoparse=True):
        self.filename = filename
        self.out_format = out_format
        self.iprlookup = iprlookup
        self.goterms = goterms
        self.autoparse = autoparse
        self.parsed = False
        self.results = list()
        self.results_byquery = dict()

        if self.autoparse:
            self.parse_results()

    def __repr__(self, ):
        if self.parsed:
            return "<InterproResultFile: {0}, {1} results for {2} queries>".format(self.filename, len(self.results), len(self.results_byquery))
        else:
            return "<InterproResultFile: {0} (UNPARSED)>".format(self.filename)
    
    def __iter__(self, ):
    # Override to iterate over results list when treated as an iterable
        return iter(self.results)

    def __getitem__(self, key):
    # Override to return all results of an ID when treated as a seekable (dict/list: IRF[seqkey])
        return self.results_byquery[key]

    def parse_results(self, ):
    # Opens and parses filename. Populates self.results and sets self.parsed=True
        if not self.results:
            self.results = list()
        if not self.results_byquery:
            self.results_byquery = dict()

        handle = open(self.filename)
        for line in handle:
            # Using list "unpacking" to auto-assign list elements to method params
            fields = line.rstrip().split('\t')
            ir = InterproResult(*fields)
            self.results.append(ir)
            try:
                self.results_byquery[ir.query_id].append(ir)
            except KeyError:
                self.results_byquery[ir.query_id] = [ir]
        self.parsed = True

    def get_query_groups(self, ):
    #DEPRECATED - use dictionary/indexable treatment or self.results_byquery'''
        pass
    # Returns a list of lists, where the element lists are all InterproResult objects with same query_id
        if not self.parsed:
            print "Parse results file before attempting grouping"
            return None
        groups = list()
        tmp_group = list()
        tmp_group.append(self.results[0])
        for result in self.results[1:]:
            if result.query_id == tmp_group[0].query_id:
                tmp_group.append(result)
            else:
                groups.append(tmp_group)
                tmp_group = list()
                tmp_group.append(result)
        groups.append(tmp_group)
        return groups


class InterproResult():
    def __init__(self, 
                 query_id, 
                 checksum, 
                 query_length, 
                 method, 
                 result, 
                 result_desc, 
                 match_start, 
                 match_end, 
                 evalue, 
                 status, 
                 run_date, 
                 ipr=None, 
                 ipr_desc=None,
                 go_desc=None):
        self.query_id = query_id
        self.checksum = checksum
        self.query_length = query_length
        self.method = method
        self.result = result
        self.result_desc = result_desc
        self.match_start = match_start
        self.match_end = match_end
        self.evalue = evalue
        self.status = status
        self.run_date = run_date
        self.ipr = ipr
        self.ipr_desc = ipr_desc
        self.go_desc = go_desc
        self.go_terms = self._parse_go_desc()
    
    def _parse_go_desc(self, ):
    # Returns a list of GO:NN..N terms as strings, or None
        import re
        if self.go_desc == None:
            return None
        go_pat = r"(?P<go_term>GO:[0-9]+)"
        terms = re.findall(go_pat, self.go_desc)
        if terms == []:
            return None
        return terms

    def __repr__(self, ):
        return "<InterproResult: {0}: {1}\t{2}\t{3} ({4})>".format(self.query_id, self.method, self.result, self.ipr, self.ipr_desc)
    
    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self.__eq__(other)
     

class Interpro():
    
    def __init__(self, sequence_id, results_dir='~/hpf_interpro', appl_str=None, seqtype='p', auto_create=True, cleanup=True):
        self.sequence_id = sequence_id
        self.results_dir = os.path.expanduser(results_dir)
        self.appl_str    = appl_str
        self.seqtype     = seqtype
        self.cleanup     = cleanup

        self.sequence_file = os.path.join(self.results_dir, str(self.sequence_id) + ".fasta")
        self.outfile = os.path.join(self.results_dir, str(self.sequence_id) + ".interpro.out")

        if auto_create:
            self._check_dirs()
            self._make_seqfile()
    
    def _check_exec(self, ):
        command = ['which', IPRSCAN]
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError as e:
            print "Interpro executable '{0}' does not exist on Path".format(IPRSCAN)
            raise

    def _check_dirs(self, ):
        if not ( os.path.exists(self.results_dir) and os.path.isdir(self.results_dir) ):
            raise OSError("Results directory {0} does not exist. Create and re-run".format(self.results_dir))

    def _make_seqfile(self, ):
        from hpf.hddb.db import Session, Sequence
        session=Session()
        sequence = session.query(Sequence).get(self.sequence_id)
        if not sequence:
            raise Exception("Getting sequence object from database failed")
        outhandle = open(self.sequence_file, 'w')
        outhandle.write(">hpf_seqid|{0}\n".format(sequence.id))
        outhandle.write("{0}\n".format(sequence.sequence))
        outhandle.close()
        session.close()

    def run(self, ):
        self._check_exec()
        command = [IPRSCAN, "-cli", "-i", self.sequence_file, "-o", self.outfile, "-format", "raw", "-seqtype", self.seqtype, "-verbose", "-iprlookup", "-goterms"]
        if self.appl_str:
            command.append("-appl")
            command.append(self.appl_str)
        print "Running interpro command: {0}".format(command)
        try:
            subprocess.check_call(command)
        except OSError as e:
            print "Execution failed, raising exception:\n"
            raise(e)
        if self.cleanup:
            self._cleanup()

    def _cleanup(self, ):
        # Remove sequence file
        try:
            os.unlink(self.sequence_file)
        except:
            print "Removing tmp file '{0}' failed. Leaving file".format(self.sequence_file)
        

