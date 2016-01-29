import os
import subprocess
import re
import decimal
import tempfile
from Bio.PDB.Structure import Structure
from hpf.runtime import debug

class Mammoth(object):

    def __init__(self,cl,parse=True):
        """
        @param parse: If false, return None after running. 
        """
        self.cl = cl
        self.parse = parse
        
    def run(self):
        return do_alignment(self.cl,parse=self.parse)

def do_alignment(cl,parse=True):
    command = str(cl)
    if cl.tempd == None:
        cwd = tempfile.mkdtemp()
    else:
        cwd = cl.cwd
    
    debug("Running mammoth",command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=open(os.devnull,"w"),cwd=cwd)
    if not parse:
        return p.communicate()
    
    # Now parse output
    if cl.output:
        # Consume the process' pipe so we don't block
        for l in p.stdout:
            pass
        parse_output = open(cl.output)
    else:
        parse_output = p.stdout
    ini_psi, ini_rms, end_psi, end_rms, zscore, evalue, pred_seq, pred_ss, exp_seq, exp_ss = _parse_mammoth(parse_output)
    # Wait for process to finish
    rcode = p.wait()
    if rcode != 0:
        raise subprocess.CalledProcessError(command,rcode)
    
    mm = MammothAlignment(ini_rms=ini_rms, ini_psi=ini_psi, end_rms=end_rms, end_psi=end_psi, zscore=zscore,
        evalue=evalue, pred_seq=pred_seq, pred_ss=pred_ss, exp_seq=exp_seq, exp_ss=exp_ss)
    mm._cwd = cwd
    return mm
    
class MammothAlignment:
    """Represents a Mammoth Alignment.
    Has variables ini_rms, ini_psi, end_rms, end_psi"""
    
    def __init__(self, ini_rms=None, ini_psi=None, end_rms=None, end_psi=None, zscore=None, evalue=None,
            pred_seq=None, exp_seq=None, pred_ss=None, exp_ss=None):
        self.ini_rms = ini_rms
        self.ini_psi = ini_psi
        self.end_rms = end_rms
        self.end_psi = end_psi
        self.zscore = zscore
        self.evalue = evalue

        #kdrew: adding in final structural alignment strings
        self.pred_seq = pred_seq
        self.exp_seq = exp_seq
        self.pred_ss = pred_ss
        self.exp_ss = exp_ss

    def __str__(self):
        #+"\npred_seq: "+self.pred_seq+"\npred_ss:  "+self.pred_ss+"\nexp_ss:   "+self.exp_ss+"\nexp_seq:  "+self.exp_seq
        return "ini_rms:"+str(self.ini_rms)+" ini_psi:"+str(self.ini_psi)+" end_rms:"+str(self.end_rms)+" end_psi:"+str(self.end_psi) 
        
    def none(self):
        return (self.ini_psi==None) or (self.ini_rms==None) or (self.end_psi==None) or (self.end_rms==None)
    
    def zipped(self):
        """
        Zip the sequence indices corresponding to the structural alignment
        @return: Tuple list of (predicted_index,experiment_index), 
            ie. [(None,0),(0,1)....]
        """
        def index(seq_str):
            x=0
            for res in seq_str:
                if res != '.':
                    yield x
                    x+=1
                else:
                    yield None
        
        return zip(index(self.pred_seq),index(self.exp_seq))
        
            

class MammothCL:
    """Mammoth Command Line parameters"""
    #"mammoth -e %s -p %s -r 1 -o mammoth.log"
    
    def __init__(self, experiment,prediction, command="mammoth",cwd=None,tempd=True,output=None,level=1,verbose=1):
        self.command = command
        self.cwd = cwd
        self.tempd = tempd
        self._temp_files = []
        self.experiment = self._temp(experiment)
        self.prediction = self._temp(prediction)
        self.output = output
        self.level = level
        self.verbose = verbose

    
    def set_output(self, output_file, level=None, verbose=None):
        self.output = output_file
        if level:
            self.level = int(level)
        if verbose:
            self.verbose = int(level)

    def _temp(self,structure):
        if isinstance(structure,basestring):
            assert(os.path.exists(structure))
            return structure
        elif isinstance(structure,Structure):
            temp = tempfile.NamedTemporaryFile("w")
            # Save the structure in a tempfile
            from Bio.PDB.PDBIO import PDBIO
            io = PDBIO()
            io.set_structure(structure)
            io.save(temp)
            self._temp_files.append(temp)
            return temp.name
        else:
            print type(structure)
            raise TypeError(structure)

    def __del__(self):
        for temp in self._temp_files:
            try:
                temp.close()
            except:
                pass
    
    def __str__(self):
        cline = self.command
        cline += " -e %s -p %s" % (self.experiment,self.prediction)
        if self.output:
            cline += " -o %s" % self.output
        if self.level != None:
            cline += " -r %i" % self.level
        if self.verbose != None:
            cline += " -v %i" % self.verbose
        return cline


    
def _parse_mammoth(open_file):
    """Parse a mammoth output file"""
    rms_regex = re.compile("RMS= *[0-9]+.[0-9]+")
    psi_regex = re.compile("PSI\(.*\)= *[0-9]+.[0-9]+")
    def value(regex, line):
        match = regex.search(line)
        return line[match.start():match.end()].split("=")[-1].strip()
    def psi_rmsd(line):
        score = value(rms_regex, line)
        psi = value(psi_regex, line)
        return decimal.Decimal(psi), decimal.Decimal(score)

    # PSI(ini)=   98.75  NALI=  79  NORM=  80  RMS=    4.45  NSS=  60
    # PSI(end)=   92.50  NALI=  74  NORM=  80  RMS=    3.78
    # Sstr(LG)=  890.85  NALI=  74  NORM=  80  RMS=    3.78
    #E-value=     0.17402446E-04
    #Z-score=      11.254682        -ln(E)=      10.958900

    try:
        ini_psi, ini_rmsd, end_psi, end_rmsd, zscore, evalue = (None,None,None,None,None,None)
        pred_seq = ""
        pred_ss = ""
        exp_seq = ""
        exp_ss = ""
        for line in open_file:
            if re.compile(" PSI\(ini\).*").search(line):
                ini_psi, ini_rmsd = psi_rmsd(line)
            elif re.compile(" PSI\(end\).*").search(line):
                end_psi, end_rmsd = psi_rmsd(line)
            elif re.compile("Z\-score\=.*").search(line):
                zscore = decimal.Decimal(line.split()[1])
            elif re.compile("E\-value\=.*").search(line):
                evalue = decimal.Decimal(line.split()[1])
            #kdrew: add in final structural alignment parse
            elif re.compile("^Prediction ").search(line):
                #kdrew: order of output from mammoth
                if len(pred_seq) == len(pred_ss):
                    pred_seq+=line.partition("Prediction ")[2].replace(' ','').strip()
                else:
                    pred_ss+=line.partition("Prediction ")[2].replace(' ','').strip()
            elif re.compile("^Experiment").search(line):
                #kdrew: order of output from mammoth
                if len(exp_seq) == len(exp_ss):
                    exp_ss+=line.partition("Experiment ")[2].replace(' ','').strip()
                else:
                    exp_seq+=line.partition("Experiment ")[2].replace(' ','').strip()

        return ini_psi, ini_rmsd, end_psi, end_rmsd, zscore, evalue, pred_seq, pred_ss, exp_seq, exp_ss
    finally:
        open_file.close()
    return None, None, None, None, None, None, None, None, None, None


