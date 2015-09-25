#!/usr/bin/env python

# Defines DenovoResultFile and DenovoResult classes, made to hold results
# returned from Rosetta denovo on the WCGrid. Files are in "silent" format
# (many decoy records per file).

import re

class DenovoResultFile():
    """Defines a class for representing a Rosetta de Novo structure prediction result file
    (a silent-format file with header info and many silent records representing decoys). Stores
    inidividual silent records as a list and dict (by index) of DenovoResult objects. Defines
    iterable and hashable treatment to access results list and dict, respectively.
    NOTE that the parse() method populates instance variable and does not return results.
    NOTE that the denovo results are 1-INDEXED (to fit the rosetta clusterer index scheme)
    """
    
    def __init__(self, filename, prediction_code, autoparse=True):
        """
        Note: a num_results value of -1 indicates that the results file has not yet been parsed
        num_results value of 0 indicates that there are no results in the given file
        """
        self.filename = filename
        self.prediction_code = prediction_code
        self.num_results = -1
        self.discarded = 0
        self.sorted = False
        self.results = list()
        self.results_dict = dict()

        self.sequence, self.score_key = self._parse_header()
        self.sequence_length = len(self.sequence)

        if autoparse: self.parse()
            
    def __repr__(self, ):
        return "<DenovoResultFile: {0}, {1} (results: {2})>".format(self.filename, self.prediction_code, self.num_results)
   
    def __iter__(self, ):
        """Treated as an iterable, iterates over the results list"""
        return iter(self.results)

    def __getitem__(self, key):
        """
        Treated as a hashable, returns DenovoResult obj from results_dict (key is DenovoResult.index, as id is NOT UNIQUE)
        NOTE: THE DenovoResult.id FROM ROSETTA SILENT FILES IS NOT UNIQUE
        """
        return self.results_dict[key]

    def _parse_header(self, ):
        """Parses the sequence and score key line from the given file. None returned if not found"""
        handle = open(self.filename)
        sequence, score_key = None, None
        for line in handle:
            fields = line.rstrip().split()
            if fields[0] == 'SEQUENCE:':
                sequence = fields[1]
            elif fields[0] == 'SCORE:' and fields[1] == 'score':
                score_key = line.rstrip()
            if sequence and score_key:
                break
        handle.close()
        return sequence, score_key

    def parse(self, ):
        """Populates self.results and self.results_dict, a list and dict of DenovoResult objects"""
        handle = open(self.filename)
        index = 1
        id, score, score_line = None, None, None
        lines = list()
        for line in handle:
            if re.search(r"SCORE:\s+-?[0-9]+\.[0-9]+|SCORE:\s+nan", line):
                # Create DenovoResult for the previous record and add it to results (if there is a prev. record)
                if id and score != None and lines != []:
                    if self._check_record(lines):
                        self._add_record(id, index, score, score_line, lines)
                        index += 1
                    else:
                        print "Record {0} corrupt or incomplete. Discarding".format(id)
                        self.discarded += 1
                # Set values for the new record
                id = line.rstrip().split()[-1]
                score = line.rstrip().split()[1]
                score_line = line.rstrip()
                if re.match(r"-?[0-9]+\.[0-9]+", score):
                    score = float(score)
                else:
                    score = 0.0
                lines = [line]
            elif re.search(r"\s+[0-9]+\s+[a-zA-Z]\s+-?[0-9]+\.[0-9]+", line):
                lines.append(line)
            elif re.match(r"SEQUENCE:\s+[A-Z]+", line):
                continue
            elif re.match(r"SCORE:\s+score", line):
                continue
            elif line.rstrip() == "":
                continue
            elif re.match(r"Structures Completed:", line):
                continue
            elif re.match(r"Total Attempts:", line):
                continue
            elif re.match(r"Attempts Per Structure:", line):
                continue
            else:
                #raise Exception("Line '{0}' not recognized".format(line.rstrip()))
                print "Line '{0}' not recognized. Skipping..".format(line.rstrip())
                continue
        # Add last record
        if self._check_record(lines):
            self._add_record(id, index, score, score_line, lines)
        else:
            print "Record {0} corrupt or incomplete. Discarding".format(id)
            self.discarded += 1
        self.num_results = index
        handle.close()
        print "Parsing {0} completed successfully. {1} results after {2} discarded".format(self.filename, self.num_results, self.discarded)
    
    def _add_record(self, id, index, score, score_line, lines):
        """A function to create and add a DenovoResult record obj to the instance (results list and dict)"""
        result = DenovoResult(id, index, score, score_terms=score_line, lines=lines)
        self.results.append(result)
        self.results_dict[index] = result

    def _check_record(self, lines):
        """A function to check the integrity of a DenovoResult silent record before
        adding it to the instance's results list and dict. Currently, simply checks 
        that the number of residue lines (lines - 1 to subtract score line) in the 
        parsed lines list is equal to the length of the sequence (must be to be valid).
        Can be re-implement for a more sophisticated check..
        """
        if len(lines) - 1 == self.sequence_length:
            return True
        return False

    def sort_by_score(self, ):
    # Sorts results list by score, lowest (best) to highest (worst).
    # Sorts inplace, does not return list
        if self.num_results < 0:
            raise Exception("Results for file '{0}' have not yet been parsed".format(self.filename))
        self.results.sort(key=lambda result: result.score)
        self.sorted = True

    def get_best(self, ):
    # Returns DenovoResult object with best (lowest) score
        if self.num_results < 0:
            raise Exception("Results for file '{0}' have not yet been parsed".format(self.filename))
        if not self.sorted:
            self.sort_by_score()
        return self.results[0]
    
    def get_worst(self, ):
    # Returns DenovoResult with worst (highest) score
        if self.num_results < 0:
            raise Exception("Results for file '{0}' have not yet been parsed".format(self.filename))
        if not self.sorted:
            self.sort_by_score()
        return self.results[-1]
    
    def get_top_count(self, count):
    # Returns the best 'count' (number) of results as a list of DenovoResult objects
        if count > self.num_results:
            print "Requested more results than available. Returning all results"
            count = self.num_results
        if not self.sorted:
            self.sort_by_score()
        return self.results[:count]

    def get_top_percent(self, percent):
    # Returns the best 'percent' percent of results as a list of DenovoResult objects
        if percent > 100 or percent < 0 :
            raise Exception("Percent out of range 0-100. Seriously.")
        if not self.sorted:
            self.sort_by_score()
        return self.results[:int(self.num_results * (percent / 100.0))]

    def print_results(self, ):
        print "{0}, code: {1}".format(self.filename, self.prediction_code)
        print "SEQUENCE: {0}".format(self.sequence)
        print self.score_key
        if self.num_results < 0:
            print "Results have not yet been parsed"
        elif self.num_results == 0:
            print "No decoy results in given file"
        for result in self.results:
            result.print_result()

    def write_to_file(self, outfile=None, percent=None, count=None, results_list=None):
    # Writes a silent file of denovo result records. Can do by top percent, count, or a given list of DenovoResult objs
        if not outfile:
            outfile = "{0}.denovo.result".format(self.prediction_code)
        outhandle = open(outfile, 'w')
        if percent:
            results = self.get_top_percent(percent)
        elif count:
            results = self.get_top_count(count)
        elif results_list:
            results = results_list
        else:
            results = self.results
        outhandle.write("SEQUENCE: {0}\n".format(self.sequence))
        outhandle.write("{0}\n".format(self.score_key))
        for result in results:
            result.write_to_filehandle(outhandle)
        outhandle.close()
        return outfile

    def add_silent_result(self, ):
        raise Exception("Not yet implemented")
        self.sorted = False

    def remove_result(self, id):
        raise Exception("Not yet implemented")


class DenovoResult():
    """Represents a single deNovo result silent record. Parameters:
    id    - Rosetta decoy ID (last column of all entry lines)
    index - index of decoy in file. NOTE: STARTS AT 1 (not 0)
    score - Decoy score (first value-column of SCORE line)
    score_terms   - The decoy's 'SCORE: <values>' line 
    lines - List of all decoy's lines (including SCORE line)
    """
    def __init__(self, id, index, score, score_terms=None, lines=[]):
        self.id = id
        self.index = index
        self.score = score
        self.lines = lines
        if score_terms:
            self._set_score_terms(score_terms)

    def __repr__(self, ):
        return "<DenovoResult ID: {0} (index {1}), score: {2}>".format(self.id, self.index, self.score)
    
    def _set_score_terms(self, score_terms):
        """Parses score terms from score_terms line according to this Rosetta score key: 
           SCORE: score env pair vdw hs ss sheet cb rsigma hb_srbb hb_lrbb rg co contact rama bk_tot  \
           fa_atr fa_rep fa_sol h2o_sol hbsc fa_dun fa_intra fa_pair fa_plane fa_prob fa_h2o \ 
           h2o_hb gsolt sasa omega_sc description
        Currently only sets a few instance vars (total: 32 fields)
        """
        terms = score_terms.split()
        self.radius_gyration = float(terms[12])
        self.contact_order   = float(terms[13])
        self.contact         = float(terms[14])
        self.rama            = float(terms[15])

    def print_result(self, ):
        for line in self.lines:
            print line.rstrip()
    
    def write_to_filehandle(self, outhandle):
        for line in self.lines:
            outhandle.write(line)

    def add_lines(self, lines):
        """Can take either a list or a single string. Adds lines to end of record"""
        self.lines += lines
        
