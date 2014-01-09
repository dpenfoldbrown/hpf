'''
Created on Oct 12, 2009

@author: Patrick
'''
import os
import subprocess
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from hpf.runtime import runtime

def percent_alpha_beta(sequence, ss):
    """A function to calculate the percent alpha and beta of a sequence given the sequence and 
    Psipred SS prediction string (eg 'CCEECCCCCEEEEECCCCCCCCCC')
      sequence    - (str) the sequence in question
      ss          - (str) a string representing the secondary structure prediction of the given sequence
    Returns:
      sequence length, percent_aplha, percent_beta
    """
    if len(sequence) != len(ss):
        raise Exception("Sequence ({0}) does not match SS pred ({1})".format(sequence, ss))
    ss = ss.upper()
    length = len(sequence)
    percent_alpha = float(ss.count('H')) / float(length)
    percent_beta  = float(ss.count('E')) / float(length)
    return length, percent_alpha, percent_beta


class PsipredAlphabet(SingleLetterAlphabet):
    # strand, helix, loop
    letters = "EHC"


class PsipredPrediction(Seq):
    
    def __init__(self, pred, weights):
        """
        @type pred: str. Is a string of the SS prediction (eg: CCCHHEEC)
        @type weights: list<int>
        """
        Seq.__init__(self, pred, PsipredAlphabet())
        self.prediction = pred
        self.weights = weights
        
    def label(self, num):
        return PsipredAlphabet.letters.index(self[num])


class HPFPsipredWrap():
    """
    A complete wrapper for running Psipred for HPF applications and storing results
    in HPF DB (hpf.sequencePsiPred). Takes sequence key, will create fasta, checkpoint,
    and psipred files, and will push psipred result by sequence key and ginzu_version 
    to DB.
    Note: In supplying a ginzu version, choose the ginzu version that contains the
    version of Psipred that you run (eg, Psipred 3.2 is ginzu_version 4)
    The default value for ginzu_version should reflect the version of psipred of the 
    system on which this code is run.
    Note: Requires the location of the NR blast database to make checkpoint file
    """
    def __init__(self, sequence_key, nr_db, ginzu_version="4", dir=None, autorun=True, dbstore=True, debug=True):
        """Variables:
        self.prediction - The PsipredPrediction object returned from Psipred32(...).run()
        self.dbo        - If dbstore=True, the Psipred ORM object (DataBaseObject) from the HPF DB
        """
        from hpf.hddb.db import Session, Sequence

        self.sequence_key  = sequence_key
        self.nr_db = os.path.abspath(os.path.expanduser(nr_db))
        self.ginzu_version = ginzu_version
        self.dir = dir if dir else os.getcwd()
        self.dbstore = dbstore
        self.debug = debug

	#kdrew: commenting out because "nr" is not a file but a location
        #if not os.path.isfile(self.nr_db):
        #    raise Exception("Must provide a valid NR database file")
        
        self.session  = Session()
        self.sequence = self.session.query(Sequence).get(self.sequence_key)
        if not self.sequence:
            raise Exception("No sequence with key {0} exists in DB".format(self.sequence_key))
        
        self.fasta_file = "{0}.fasta".format(self.sequence_key)
        self.chkpt_file = "{0}.chk".format(self.sequence_key)
        self.psipred_file = "{0}.psipred".format(self.sequence_key)

        # Set in run (dbo optionally)
        self.prediction = None
        self.dbo = None

        if autorun:
            self.prediction = self.run()

    def get_prediction_string(self, ):
        """Returns the prediction string from the PsipredPrediction object (if populated)"""
        if not self.prediction:
            return None
        return self.prediction.prediction

    def run(self, ):
        """
        Cobbles together elements for running Psipred (writes fasta, runs blast to get checkpoint,
        runs psipred and psipass2, and then parses results and uploads to DB if set to).
        IF dbstore is true, adds Psipred ORM object to hpf database and sets self.dbo.
        RETURNS hpf.pdb.psipred.PsipredPrediction object
        """
        import subprocess
        from Bio import SeqIO
        from Bio.Blast.NCBIStandalone import blastpgp
        from hpf.hddb.db import Psipred as PsipredORM, PsipredFactory
        from sqlalchemy.exc import IntegrityError

        # Write fasta file
        if self.debug: print "Psipred: writing fasta file"
        with open(self.fasta_file, 'w') as handle:
            SeqIO.write([self.sequence.record], handle, "fasta")

        # Get exe path and run blast
        if self.debug: print "Psipred: running blastpgp against '{0}' DB to create checkpoint file".format(self.nr_db)
        blast_cmd = subprocess.Popen(["which", "blastpgp"], stdout=subprocess.PIPE).communicate()[0].strip()
        result,error = blastpgp(blastcmd=blast_cmd, 
                                program='blastpgp', 
                                database=self.nr_db, 
                                infile=self.fasta_file, 
                                npasses=3, 
                                checkpoint_outfile=self.chkpt_file, 
                                expectation=1e-4, 
                                model_threshold=1e-4, 
                                align_outfile="/dev/null")
        # Note: must call something that blocks on blastpgp results (need to wait for cmd to finish)
        res = result.read()
        err = error.read()
        if self.debug:
            print "Result: ", res
            print "Error/Warning: ", err

        # Create Psipred options object and run Psipred 3.2 on them
        if self.debug: print "Psipred: running Psipred 3.2"
        options = PsipredOptions(self.fasta_file,
                                 profile=self.chkpt_file,
                                 output=self.psipred_file+".1",
                                 output2=self.psipred_file+".2",
                                 horiz=self.psipred_file,
                                 cwd=self.dir)
        self.prediction = Psipred32(options).run()
        
        # Add to database (optional) and return
        if self.dbstore:
            if self.debug: print "Psipred: adding psipred DBO to hpf database"
            psipred_dbo = PsipredFactory().create(self.prediction, sequence_key=self.sequence.id, ginzu_version=self.ginzu_version)
            self.session.add(psipred_dbo)
            try:
                self.session.commit()
            except IntegrityError:
                print "Psipred entry for <seq {0}, ginzu_version {1}> already exists in DB. Returning existing object".format(self.sequence_key, self.ginzu_version)
                self.session.rollback()
                psipred_dbo = self.session.query(PsipredORM).filter_by(sequence_key=self.sequence.id, ginzu_version=self.ginzu_version).first()
            self.session.refresh(psipred_dbo)
            self.dbo = psipred_dbo
        return self.prediction


class PsipredWriter(object):
    """
    Writes a psipred prediction file.
    This doesn't implement the multiples of 10 column numbers.
    """


#Conf: 973100001489875201002235666777640797544278888877888878743223
#Pred: CCCCHHHHCCCCCCCCEEEHHHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#  AA: MFRRAAMAAAAPSDSEVRVQKVDKLDLVFNILTKPPVYGAGKGNNPPKAPAPRRPAATGG
#              10        20        30        40        50        60
#
#
#Conf: 66666676777788544143311330247889999843789
#Pred: CCCCCCCCCCCCCCCEEEHHHHCCCHHHHHHHHHHHHCCCC
#  AA: DHGSGGAVAGRKQPGVVSIEDINKRSENYIRDRKRMFFGQN
#              70        80        90       100
    
    def write(self, handle, prediction, record):
        conf = prediction.weights
        pred = str(prediction)
        aa = str(record.seq)
        length = len(pred)/60
        if length*60 < len(pred):
            length+=1
        def _line(i):
            start = (i)*60
            end = start+60
            
            print >>handle, ""
            print >>handle, "Conf: "+"".join([str(w) for w in conf[start:end]])
            print >>handle, "Pred: "+pred[start:end]
            print >>handle, "  AA: "+aa[start:end]
            print >>handle, ""
            
        print >>handle, "# PSIPRED HFORMAT (PSIPRED V2.6 by David Jones)"
        for i in range(length):
            _line(i)
    

class PsipredOptions(object):

    def __init__(self, fasta, profile=None, output=None, output2=None, horiz=None, cwd=None):
        """
        @param fasta: sequence FASTA file
        @param profile: BLAST checkpoint profile
        """
        self._temp = None
        self.fasta = fasta
        self.profile = profile
        self.single = (profile == None)
        self.output = output
        self.output2 = output2
        self.horiz = horiz
        self.cwd = cwd if cwd else os.getcwd()
        
        assert os.path.exists(fasta)
        if not self.single:
            assert os.path.exists(profile)

    def _mtx(self, name=None):
        """
        echo $tmproot.chk > $tmproot.pn
        echo $tmproot.fasta > $tmproot.sn
        $ncbidir/makemat -P $tmproot
        
        or
        
        $execdir/seq2mtx $1 > $tmproot.mtx
        """
        if not name:
            name = ".".join(os.path.basename(self.profile).split(".")[:-1])
        mtx = "%s.mtx" %name

        # Either generate the matrix from a profile or sequence alone
        if self.single:
            cmd = "seq2mtx %s > %s" % (self.fasta, mtx)
            #print cmd
            runtime().debug(cmd)
            subprocess.check_call(cmd, shell=True, cwd=self.cwd)
        else:
            for file,link in [(self.profile,name+".pn"),(self.fasta,name+".sn")]:
                subprocess.check_call("echo %s > %s" % (file,link), shell=True, cwd=self.cwd)
            makemat = subprocess.Popen("which makemat", shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
            cmd = "%s -P %s" % (makemat, name)
            #print cmd
            runtime().debug(cmd)
            subprocess.check_call(cmd,shell=True, cwd=self.cwd)
        return mtx


class Psipred32():
    """
    a run wrapper for version 3.2-style Psipred: one less weights file, single run difference
    is in how matrix is created and does not affect psipred or psipass2 commands or weights files
    used (as Psipred 3.2+ does not have 'single-run' weights files).
    See run scripts provided with psipred distribution for more information
    """
    def __init__(self, options):
        self.options = options
    
    def run(self, ):
        mtx = self.options._mtx()
        output = self.options.output
        psipred = subprocess.Popen("which psipred", stdout=subprocess.PIPE, shell=True).communicate()[0]
        # Cut the /bin/psipred off
        bin_dir = os.path.split(psipred)[0]
        root_dir = os.path.split(bin_dir)[0]
        data = os.path.join(root_dir,"data")
        
        # Single run (from runpsipred_single script):
        #psipred $tmproot.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $rootname.ss
        # Normal run (from runpsipred script):
        #psipred $tmproot.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $rootname.ss
        cmd = "psipred {0} {1}/weights.dat {1}/weights.dat2 {1}/weights.dat3 > {2}".format(mtx, data, output)
        
        runtime().debug(cmd)
        subprocess.check_call(cmd, shell=True, cwd=self.options.cwd)
        
        output2 = self.options.output2
        horiz = self.options.horiz
        
        # psipass2 command (from runpsipred script):
        #psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $rootname.ss2 $rootname.ss > $rootname.horiz
        cmd = "psipass2 {0}/weights_p2.dat 1 1.0 1.0 {1} {2} > {3}".format(data, output2, output, horiz)
        runtime().debug(cmd)
        subprocess.check_call(cmd, shell=True, cwd=self.options.cwd)
        
        # Note: Output format between Psipred versions 2.5 and 3.2 is the same (woo)
        with open(self.options.horiz) as handle:
            pred = parse(handle)
        return pred


class Psipred(object):
    
    """
    dpb: This works for old version of Psipred (-3.2)

    pn blast checkpoint file
    sn fasta file
    `makemat` NCBI toolkit
    
    $execdir/psipred $tmproot.mtx 
        $datadir/weights.dat 
        $datadir/weights.dat2 
        $datadir/weights.dat3 
        $datadir/weights.dat4 > $rootname.ss
        
    or
    
    $execdir/psipred $tmproot.mtx 
        $datadir/weights_s.dat 
        $datadir/weights_s.dat2 
        $datadir/weights_s.dat3 > $rootname.ss
        
    $execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $rootname.ss2 $rootname.ss >
        $rootname.horiz
    """
    
    def __init__(self, options):
        self.options = options
        
    def run(self, twopass=True):
        mtx = self.options._mtx()
        output = self.options.output
        psipred = subprocess.Popen("which psipred", stdout=subprocess.PIPE, shell=True).communicate()[0]
        # Cut the /bin/psipred off
        bin_dir = os.path.split(psipred)[0]
        root_dir = os.path.split(bin_dir)[0]
        data = os.path.join(root_dir,"data")
        if self.options.single:
            cmd = """psipred %s 
                    %s/weights_s.dat 
                    %s/weights_s.dat2 
                    %s/weights_s.dat3 > %s""" % (mtx,data,data,data,output)
        else:
            cmd = """psipred %s %s/weights.dat %s/weights.dat2 %s/weights.dat3 %s/weights.dat4 > %s""" % (mtx,data,data,data,data,output)
        #print cmd
        runtime().debug(cmd)
        subprocess.check_call(cmd, shell=True, cwd=self.options.cwd)
        
        output2 = self.options.output2
        horiz = self.options.horiz
        cmd = """psipass2 %s/weights_p2.dat 1 1.0 1.0 %s %s > %s""" % (data, output2, output, horiz)
        runtime().debug(cmd)
        subprocess.check_call(cmd, shell=True, cwd=self.options.cwd)
        with open(self.options.horiz) as handle:
            pred = parse(handle)
        return pred


def parse(file_handle):
    # 'seq' is the prediction sequence, not the AA sequence (eg: CCCHHEEC...)
    import re
    result_pat = r"Conf:\s[0-9]+$|Pred:\s[A-Z]+$"
    
    seq = ""
    weights = []
    
    for line in file_handle:
        if re.match(result_pat, line):
            if line.startswith("Conf"):
                weights += [int(w) for w in line.strip().split()[1]]
            elif line.startswith("Pred"):
                seq += line.strip().split()[1]
    
    assert len(seq) == len(weights)
    return PsipredPrediction(seq, weights)

