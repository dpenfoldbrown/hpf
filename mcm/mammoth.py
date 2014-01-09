'''
Necessary mammothing tools for MCM.
Created on Apr 28, 2010
@author: Patrick

@todo: Computing contact order of mammoth matches
    could and should be done here instead of relying
    on an out-of-date modified version of mammoth.
'''
import os
from hpf.pdb.mammoth import Mammoth, MammothCL
from hpf.utilities.temp import ControlFile
from hpf.runtime import runtime

PREDICTION_INDEX = "p"
EXPERIMENT_INDEX = "e"
ZSCORE = "Zscore"
LN = "-ln(E)"
EVALUE= "Evalue"
SCORE = "score"
NUM_E = "#e"
NUM_P = "#p"
NSUP = "nsup"
NSS = "nss"
PSI1 = "psi1"
PSI2 = "psi2"
CO_P = "CO_p"
CO_E = "CO_e"
PREDICTION = "prediction"
EXPERIMENT = "experiment"

#kdrew: occasionally, without obvious reason, mammoth produces '******' for a zscore, for now set equal to 0
def zeroStar(toks):
    import re
    if re.match("\*+",toks[0]):
        return 0.0

class MammothScore(object):
    """
    Defines a mammoth score as a dict-like object.  Keys are defined as
        constants by the module or the mapped names defined by the class.
    """
    
    # The constants above cannot be referenced directly because of python syntax
    # or the capitalization is weird so they are mapped to more variable 
    # friendly names upon assignment.
    variable_mapping = {PREDICTION_INDEX: "prediction_index",
                        EXPERIMENT_INDEX: "experiment_index",
                        ZSCORE: "zscore",
                        LN: "ln_e",
                        EVALUE: "evalue",
                        NUM_E: "num_e",
                        NUM_P: "num_p",
                        CO_P: "prediction_contact_order",
                        CO_E: "experiment_contact_order"
                        }
    
    def __init__(self,**kwargs):
        for key in kwargs:
            self[key] = kwargs[key]
            
    def __repr__(self, ):
        str = "<MammothScore"
        if self.__dict__.has_key('prediction_index'):
            str += " Prediction: {0},".format(self.__dict__['prediction_index'])
        if self.__dict__.has_key('experiment_index'):
            str += " Experiment: {0},".format(self.__dict__['experiment_index'])
        if self.__dict__.has_key(PREDICTION):
            str += " Pfile: {0},".format(self.__dict__[PREDICTION])
        if self.__dict__.has_key(EXPERIMENT):
            str += " Efile: {0}".format(self.__dict__[EXPERIMENT])
        str += ">"
        return str
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value
        if MammothScore.variable_mapping.has_key(key):
            alternate_key = MammothScore.variable_mapping[key]
            self.__dict__[alternate_key] = value


class MammothRun():
    """
    A Class to encapsulate (wrap) all functionality to run Mammoth. Run with custom contact order cmdline
        decoys      - list of decoy pdb filenames (decoy files must exist in directory def: CWD)
        experiment_listfile - the 'list.mammoth' file for mammoth to query predictions against
        prediction_listfile - a file to write decoy filenames (in list.mammoth format) to, to run mammoth on
        directory   - the directory where the prediction decoy files are located (default: CWD)
        outfile     - file to hold mammoth results and data
        executable  - mammoth executable
        autorun     - calls run method automatically. Set to False and run manually if desired
    Subclass to change exit checking, define custom or different commandline arguments, etc
    NOTE: run() runs mammoth program, and parse() parses the run's output and returns MammothScore objs
    """
    #TODO: Build in clean-up if necessary
    
    def __init__(self, decoys, experiment_listfile, prediction_listfile="prediction.mammoth", \
                 directory=None, outfile="mammoth.results", executable="mammoth", \
                 autorun=True, debug=False):
        self.decoys = decoys
        self.experiment_listfile = experiment_listfile
        self.prediction_listfile = prediction_listfile
        self.outfile = outfile
        self.executable = executable
        self.debug = debug
        
        self.directory = os.path.abspath(os.path.expanduser(directory)) if directory else os.getcwd()
        
        # Set when parse() is run
        self.scores = None
        
        if autorun:
            self.run()

    def _check_init(self, ):
        """Check all parameters required for Mammoth function"""
        if len(self.decoys) < 1:
            raise Exception("Decoys list is empty")
        if not os.path.isfile(self.experiment_listfile):
            raise EnvironmentError("Mammoth experiment list file '{0}' is not valid".format(self.experiment_listfile))
        if not os.path.isdir(self.directory):
            raise EnvironmentError("'{0}' is not a valid directory".format(self.directory))
        for file in self.decoys:
            if not os.path.isfile(os.path.join(self.directory, file)):
                raise EnvironmentError("Decoy file '{0}' does not exist in directory '{1}'".format(file, self.directory))

    def _check_exit(self, ):
        """Checks for the NORMAL_EXIT flag at the end of Mammoth output (custom hack via kdrew)"""
        with open(self.outfile) as handle:
            for line in handle: pass
            if line.find("NORMAL_EXIT") == -1:
                raise Exception("Mammoth did not exit normally (NORMAL_EXIT not found)")
    
    def _check_scores(self):
        """Check that number of SCORE lines in mammoth result file equals number of score objects parsed
        #kdrew: tests to see if score lines in data.mammoth (ie. lines starting with ':') match scores in data structure
        """
        import re
        with open(self.outfile) as handle:
            #if len(list(line for line in handle if re.match("^:",line))) != len(self.scores):
            #    raise Exception, "Not all scores from Mammoth results in '{0}' were read".format(self.outfile)
            #else:
            #    print "{0} scores read from mammoth outfile {1}".format(len(self.scores), self.outfile)
            
            num_score_lines = len(re.findall("^:", handle.read(), re.MULTILINE))
            print "MammothRun: Scores in file, {0}; Scores parsed, {1}".format(num_score_lines, len(self.scores))
            if num_score_lines != len(self.scores):
                raise Exception("Not all scores in Mammoth resultfile were read".format(self.scores))

    def run(self, ):
        # Run checks on required initial parameters
        self._check_init()

        # Open mammoth list context (populates mammoth list file) (self, file, directory, *structures)
        with MammothList(self.prediction_listfile, self.directory, *self.decoys) as list_file:
            
            # Create mastodon (mammoth w/contact order) command-line object and run Mammoth with it
            mastodon_cl = MastodonCL(contact_order=True, 
                                     experiment=self.experiment_listfile, 
                                     prediction=list_file, 
                                     cwd=os.getcwd(), 
                                     output=self.outfile, 
                                     command=self.executable, 
                                     level=0, verbose=0)
            #DEBUG
            if self.debug: print "Mammoth command: ", mastodon_cl
            ret = Mammoth(mastodon_cl, parse=False).run()
        
        self._check_exit()
        return ret

    def parse_scores(self, ):
        """Parses the instance's Mammoth result file. Returns a list of hpf.mcm.mammoth.MammothScore objects
        Returns None if outfile does not exist, and raises exception if number of results in scores list 
        does not match number of results in file
        """
        if not os.path.isfile(self.outfile):
            print "Could not parse Mammoth results: no result file '{0}' (try running first)".format(self.outfile)
            return None
        with open(self.outfile) as handle:
            self.scores = list(MammothSimpleParser().parse(handle, debug=self.debug))
        print "MammothRun: parsed {0} scores from Mammoth results file".format(len(self.scores))
        self._check_scores()
        return self.scores


class MammothSimpleParser(object):
    """A simple parser for Turbo Mammoth (MCM's included mammoth version), not using pyparsing
    NOTE: TO SUBCLASS, simply alter the order or contents of the FIELDS list, or the 
    SCORE_PATTERN string, to represent the different Mammoth results file.
    """

    # header :'     p     e    Zscore     -ln(E)     Evalue     score       #e    #p  nsup   nss   psi1      psi2       CO_p    CO_e'
    # example:':    2     1   1.891      2.247     0.1057      596.6       94   158    92    56   97.9       30.9        6.2     4.4  decoy_2691.pdb  100.pdb'
    
    # Define the pattern that matches a score line (and no other lines)
    SCORE_PATTERN = r"^:"
    
    # Define the ORDER-IMPORTANT fields list, including last two filenames, excluding leading ':'.
    # First tuple element is field name (the score dict key), second is Python type to cast the field's value to.
    FIELDS = [('prediction_index', int), 
              ('experiment_index', int),
              ('Zscore', float),
              ('-ln(E)', float),
              ('Evalue', float),
              ('score',  float),
              ('#e', int),
              ('#p', int),
              ('nsup', int),
              ('nss',  int),
              ('psi1', float),
              ('psi2', float),
              ('CO_p', float),
              ('CO_e', float),
              ('prediction', str),
              ('experiment', str)
             ]
    NUM_FIELDS = len(FIELDS)

    def parse(self, handle, debug=False):
        """Parses the Mammoth result file with static positioning"""
        import re
        for line in handle:
            if re.match(self.SCORE_PATTERN, line):
                score_dict = self.parse_score_line(line, debug=debug)
                #if debug: print "Score dict: ", score_dict
                yield MammothScore(**score_dict)
            else:
                if debug: print "Ignoring non-score line: '{0}...'".format(line[:25])

    def parse_score_line(self, line, debug=False):
        """Parse a score line from mammoth results file into a dict, {'field name' => value}"""
        
        # Get fields, discarding initial ':'. Sanity check the number of fields recvd
        line_fields = line.split()[1:]
        if len(line_fields) != self.NUM_FIELDS: 
            raise Exception("Mammoth score line '{0}' has an unexpected number of fields (should have {1})".format(line, self.NUM_FIELDS))

        # Add fields to score dict in loop and return score dict
        score_dict = dict()
        for i in range(self.NUM_FIELDS):
            try:
                # score_dict[field_name] = int/float/str(field_value)
                score_dict[self.FIELDS[i][0]] = self.FIELDS[i][1](line_fields[i])
            except ValueError:
                # Sometimes Mammoth sets Zscore to '******'. No idea why. Catch it and set Zscore to 0.0
                if self.FIELDS[i][0] == 'Zscore' and line_fields[i] == '******':
                    if debug: print "Mammoth score '{0}...' has Zscore = '******'. Setting Zscore to 0.0".format(line[:15])
                    score_dict[self.FIELDS[i][0]] = 0.0
                else:
                    raise
        return score_dict


class MastodonCL(MammothCL):
    """
    Defines the modified command line for Mastodon - a contact
    order computing version of Mammoth.
    """
    
    def __init__(self, contact_order, *args, **kwargs):
        MammothCL.__init__(self, *args, **kwargs)
        self.contact_order = contact_order

    def __str__(self):
        cline = MammothCL.__str__(self)
        if self.contact_order:
            cline += " -c 1"
        return cline
    
class MammothList():
    """
    Creates a list file as a context manager for running mammoth.  If there is 
    only one structure this will not create a list file, but return the 
    structure file name on enter.
    
    Mammoth does not accept single file lists and will error.

    Can be used as a context manager (via with), which will write the file automatically,
    or manually by instantiating and calling write()
    """
    
    def __init__(self, file, directory, *structures):
        self.file = file
        self.directory = directory
        self.structures = structures
        
        # Check parameters
        if not self.structures or len(self.structures) < 1:
            raise Exception("No structures provided to write to file")
        elif len(self.structures) == 1:
            raise UserWarning("Mammoth does not accept single-file lists and will fail (for unknown reasons)")

    def __enter__(self, ):
        return self.write()
    
    def __exit__(self, type, value, traceback):
        pass

    def write(self, ):
        with open(self.file, 'w') as handle:
            handle.write("MAMMOTH List\n")
            handle.write("{0}\n".format(self.directory))
            for structure in self.structures:
                handle.write("{0}\n".format(structure))
        return self.file


##
## DEPRECATED
##

class MammothMultiParser(object):
    """
    DPB: THIS IS DEPRECATED (BECAUSE IT IS WHACK). SIMPLE PARSER FOUND IN MammothSimpleParser
    
    Parses mammoth multi-run output where verbose is set to 0.
    This could and probably should be set up to parse any output based on a
    command line 'MammothCL' object.  For now it has one parsing expression
    that expects contact order which is not native to mammoth.
    """
 
    def parse(self,handle,debug=False):
        raise Warning("MammothMultiParser is DEPRECATED because it does not work. Please use MammothSimpleParser")

