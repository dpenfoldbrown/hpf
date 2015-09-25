#!/usr/bin/env python

import sys
import os
import subprocess
from collections import defaultdict
from hpf import hddb
from hpf.muscle import MuscleTask
from hpf.runtime import runtime, Flag, debug
from hpf.processing import processor
from hpf.utilities import OutputTask,consume
from Bio.Blast import NCBIStandalone, NCBIXML
from treecorr.alignment import AlignmentInfo, MultipleSequenceAlignment
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import pylab
from tempfile import NamedTemporaryFile
from numpy import mean
from hpf.cdhit import CDHit, CDHitOptions
from hpf.hddb.db import *

SYNCHRONOUS = "synchronous"
FORCE = "force"
DOMAIN_TYPES = ["psiblast","fold_recognition","pfam","msa","unassigned","user_defined"]
DOMAIN_COLORS = ["blue","green","red","magenta","cyan","yellow"]

Session = None

def thread_setup():
    global Session
    if Session==None:
        clear()
        e,b,s = setup()
        debug("Entering subprocess setup",s,Session)
        Session = s

def _blast(fasta):
    thread_setup()
    if runtime().opt(FORCE):
        debug("Forcing re-run")
    return blast(fasta,force=runtime().opt(FORCE))

def blast(fasta,db="hpf_protein",force=False):
    """
    Map families to the database.
    """
    runtime().debug(fasta)
    file = os.path.abspath(fasta)
    base = ".".join(os.path.basename(file).split(".")[:-1])
    name = base.replace("_"," ")
    dir = base
    subprocess.call("mkdir -p %s"%base, shell=True)
    runtime().debug(dir,fasta,base)
    cwd = os.getcwd()
    os.chdir(dir)
    
    try:
        raise Exception("This has been modified like crazy, don't run as is, make sure this is correct")
        #runtime().pushd(dir)
        #formatted_fasta = FormatFastaTask(file,base+".fasta").run(force=force)
        # Cluster the families representatives before blasting everything to HPF
        #cdhit_fasta = CdhitTask(formatted_fasta,formatted_fasta+".cdhit",identity=0.7,length=0.7).run()
        #alignment = MuscleTask(formatted_fasta, base+".aln", clwstrict=True).run(force=force)
        #alignment = FormatAlignmentTask(alignment,base+".alnf").run(force=force)
        #blast_xml, blast_chk = InputAlignmentBlastTask(formatted_fasta,alignment, db="hpf_protein").run(force=force)
        #blast_matches = BlastParseXMLTask(blast_xml, base+".hpf.fasta", 0.8, expect=1e-6).run(force=force)
        #blast_matches = BlastTask(cdhit_fasta,alignment,base+".hpf.fasta").run(force=force)
        #graphics = DomainGraphicsTask(base,blast_matches,base+".svg","","").run(force=force)
        
        blast_matches = base+".hpf.fasta"
        
        session = Session()
        family = session.query(Family).filter(Family.name==name).one()
        debug(family)
        with open(blast_matches) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                map = session.query(FamilySequence).filter(and_(FamilySequence.family_key==family.id, FamilySequence.sequence_key==int(record.id))).first()
                if map == None:
                    map = FamilySequence()
                    map.family_key=family.id
                    map.sequence_key = int(record.id)
                    debug(map)
                    session.add(map)
                
        session.commit()
        session.close()
            
        #runtime().popd()
    finally:
        os.chdir(cwd)

class CdhitTask(OutputTask):
    
    def __init__(self, fasta, output, identity=0.9, length=0.0):
        OutputTask.__init__(self, output)
        self.fasta = fasta
        self.length = length
        self.identity = identity
    
    def _do(self):
        CDHit(CDHitOptions(self.fasta, self.output, self.identity, self.length)).run()

class FormatFastaTask(OutputTask):
    """
    Formats the family fasta files by removing '*' and numbering the sequences.
    """
    
    def __init__(self,fasta,output):
        OutputTask.__init__(self,output)
        self.fasta = fasta
        
    def _do(self):
        with open(self.fasta) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        for i,record in enumerate(records):
            record.description = record.id
            record.id = str(i+1)
            record.seq = Seq(str(record.seq).replace("*", ""), record.seq.alphabet)
        with open(self.output,"w") as handle:
            SeqIO.write(records, handle, "fasta")
        return self.output

class BlastTask(OutputTask):
    """
    Blast each record in the fasta file individually and append new hits to 
    the output file.  
    @note: blastpgp hits segmentation faults when doing this with a large query fasta
    """

    def __init__(self,fasta,alignment,output,align_perc=0.8,expect=1e-6,db="hpf_protein",**kwargs):
        OutputTask.__init__(self,output)
        self.blast_exe = subprocess.Popen("which blastpgp", shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
        self.fasta = fasta
        self.alignment = alignment
        self.align_perc = align_perc
        self.expect = expect
        self.db = db
        self.kwargs = kwargs
        
    def _do(self):
        for i,fasta in enumerate(self.records()):
            try:
                with NamedTemporaryFile(prefix="blast",dir=os.getcwd()) as blast_out:
                    self.blast(fasta,blast_out.name)
                    self.parse(blast_out.name)
            except Exception, e:
                debug("Failed to blast/parse %i record in %s" % (i,self.fasta))
                debug("\t",e)
                continue

    def records(self):
        """
        Create temporary fasta files for each record.
        """
        with open(self.fasta) as handle:
            for record in SeqIO.parse(handle,"fasta"):
                with NamedTemporaryFile(prefix="fasta",dir=os.getcwd()) as temp_fasta:
                    SeqIO.write([record],temp_fasta,"fasta")
                    temp_fasta.flush()
                    yield temp_fasta.name
    
    def blast(self,fasta,output):
        """
        Blast the fasta, consume the output buffer, return the output filename
        """
        runtime().debug("Blasting %s with alignment %s using %s" %(fasta, self.alignment,self.blast_exe))
        r,e = NCBIStandalone.blastpgp(self.blast_exe, 
                                      self.db,
                                      fasta,
                                      align_infile=self.alignment,
                                      align_outfile=output,
                                      expectation=self.expect, 
                                      model_threshold=self.expect,
                                      npasses=3,
                                      nprocessors=1,
                                      **self.kwargs)
        consume(r)
        return output

    def parse(self,output):
        """
        Parse the Blast handle and append records to the output file
        """
        with open(output) as handle:
            runtime().debug("Parsing blast results")            
            blast_records = NCBIXML.parse(handle)
            matches = set()
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        n = float(min(hsp.align_length,blast_record.query_letters))
                        m = float(max(hsp.align_length,blast_record.query_letters))
                        if n/m < self.align_perc:
                            continue
                        if hsp.expect > self.expect:
                            continue
                        matches.add(alignment.title.split()[0][4:])
        runtime().debug("Found %i hits" % len(matches))
        if os.path.exists(self.output):
            with open(self.output) as fasta_file:
                already_found = set([r.id for r in SeqIO.parse(fasta_file, "fasta")])
            debug("Found %i  already" % len(matches.intersection(already_found)))
            matches = matches-already_found
        
        debug("Found %i new matches"%len(matches))
        if len(matches) > 0:
            with hddb.connect("ddbCommon") as cursor:
                records = list(hddb.proteins(cursor, sequence_key=matches))
            runtime().debug("Formatting blast results into fasta file")
            with open(self.output,"a") as fasta_file:
                SeqIO.write(records, fasta_file, "fasta")

def domain_types(records):
    """
    Split up and list relative coverage by domain_type
    """
    domain_type = defaultdict(lambda: list())
    for record in records:
        for type in record.annotations['domain_type']:
            relative_coverage = float(record.annotations['domain_type'][type])/float(len(record.seq))
            domain_type[type].append(relative_coverage)
    return domain_type

def annotate_records(fasta):
    """Attach length of domain_type coverage over each protein in fasta file"""
    with open(fasta) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        protein = {}
        for r in records:
            protein[int(r.id)] = r
            r.annotations['domain_type'] = defaultdict(lambda: 0)
        
    with hddb.connect("ddbCommon") as cursor:
        query = """select d.parent_sequence_key, d.domain_type, sum(length(sequence)) 
        from domain d join sequence s on d.domain_sequence_key=s.id 
        where parent_sequence_key in (%s) 
        group by domain_type""" % ",".join([str(id) for id in protein.keys()])
        debug(query)
        cursor.execute(query)
        for protein_key, domain_type, dlen in cursor.fetchall():
            debug(len(protein[protein_key].seq),domain_type,dlen)
            protein[protein_key].annotations['domain_type'][domain_type] = dlen
    return protein.values()

class DomainGraphicsTask(OutputTask):
    """
    Turn a families domain coverage information into a set of graphics.
    """
    
    def __init__(self, family_name, fasta, piechart, psi_hist, fr_hist, **kwargs):
        OutputTask.__init__(self, piechart, psi_hist, fr_hist)
        self.family_name = family_name
        self.fasta = fasta
        self.kwargs = kwargs
    
    def _do(self):
        records = annotate_records(self.fasta)
        domain_type = domain_types(records)
        self.piechart(domain_type)
        
    def piechart(self, domain_type):
        pie_file = self.output[0]
        pylab.figure(figsize=(5,5))
        
        probs,labels = ([],[])
        for type in DOMAIN_TYPES:
            mean_coverage = mean(domain_type[type]) if len(domain_type[type])>0 else 0.0
            debug(self.family_name,type,mean_coverage)
            probs.append(mean_coverage)
            labels.append(type)#+":%0.2f" % mean_coverage)
        
        #ax = pylab.axes([0.6, 0.6, 0.4, 0.4])
        explode = [0.05 for i in xrange(len(DOMAIN_TYPES))]
        patches, texts = pylab.pie(probs,explode=None,labels=None,shadow=False,colors=DOMAIN_COLORS)
        pylab.figlegend(patches, labels, "lower left", fancybox=True, markerscale=0.2)
        pylab.title(self.family_name)
        pylab.savefig(pie_file)
        
    
class FormatAlignmentTask(OutputTask):
    """
    Reformats a clustal alignment into a Blast PSI-Blast alignment.
    PSI-Blast has some psuedo-clustal format that removes a bunch of stuff.
    """
    
    def __init__(self, alignment, output, **kwargs):
        self.alignment = alignment
        self.output = output
        
    def _do(self):
        with open(self.alignment) as handle:
            records = len(list(AlignIO.read(handle, "clustal")))
        with open(self.alignment) as handle:
            with open(self.output,"w") as output:
                lines = handle.readlines()[2:]
                for i,line in enumerate(lines):
                    mod = (i+1)%(records+2)
                    #runtime().debug(i,mod,line[:-1])
                    if mod==records+1:
                        #print i,line[:-1]
                        continue
                    output.write(line)
            

def main(*args):
    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
    runtime().debug("Using processor",pool)
    pool.make_tasks(lambda: args)
    consume(pool.run(_blast))
    
def tasks(*args):
    if len(args) == 0:
        assert False
        #tasks = [os.path.join(dir,dir+".fas") for dir in os.listdir(os.getcwd())]
    else:
        tasks = args
    runtime().debug("TASKS",tasks)
    return tasks
    
def _do(argv):
    r = runtime()
    r.description("""
    families.py [-options] args
    Blast protein famlies against HPF database and gather representatives.
    """)
    r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    r.add_option(Flag(FORCE, "f", description="Force re-running of all tasks."))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
