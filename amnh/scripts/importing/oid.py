#!/usr/bin/env python

import sys, os
import getopt
from hpf.runtime import runtime, Flag
from hpf.utilities import consume
from hpf.processing import processor, SYNCHRONOUS
from hpf.amnh.codeml import PositiveSelectionParser
from hpf.amnh.oid import oid_hpf
from StringIO import StringIO
from hpf.amnh.align import CulledColumnMapper

oid_key = 1
bound=False

class LRTParser(object):
    def __init__(self, alignment, culled, codeml):
        self.alignment = alignment
        self.culled = culled
        self.codeml = codeml
        self.mapper = CulledColumnMapper(alignment,culled)

    def parse(self,handle):
        start = False
        for line in handle:
            if line.strip().startswith("Positively Selected Sites : Model 8"):
                start = True
                continue
            if not start:
                continue
            if line.startswith("Common"):
                break
            pieces = line.strip().split()
            if len(pieces)!=6 or not pieces[0].isdigit():
                continue
            runtime().debug("codeml",line)
            column, aa, pr, pm, _pm_, se = pieces
            # Indices start at 1 in CodeML, convert to 0
            column = column-1
            assert column >= 0
            pr = pr.split("*")[0]
            original_column = self.mapper[column]
            from hpf.hddb.db import PositiveSelection
            ps = PositiveSelection(codeml_key = self.codeml.id,
                                   column = original_column,
                                   probability = pr,
                                   post_mean = pm,
                                   stderr = se)
            yield ps
            
        
class OIDImporter(object):
    """
    Import a set of OID files into the database
    """
    def __init__(self,
                 familyName, 
                 alignFile,
                 alignColcullLog,
                 alignSeqcullLog,
                 treeFile,
                 treeDiagCharsFile,
                 codemlFile=None, 
                 alignFormat="fasta",
                 oid_key=None):
        self.familyName = familyName
        self.treeFile = treeFile
        self.treeDiagCharsFile = treeDiagCharsFile
        self.alignFile = alignFile
        self.alignColcullLog = alignColcullLog
        self.alignSeqcullLog = alignSeqcullLog
        self.codemlFile = codemlFile
        self.alignFormat = alignFormat
        self.oid_key = oid_key
        
    def merge(self):
        from hpf.hddb.db import Session,Family
        self.session = Session()
        
        self.family = self.session.query(Family).filter(Family.name==self.familyName).first()
        if not self.family:
            runtime().debug("Creating family",self.familyName)
            self._family()
            self._alignment()
            self._tree()
        else:
            self.alignment = self.family.alignment
            self.tree = self.alignment.tree
            runtime().debug("Found family",self.family.id)

        if not self.family.alignments[0].tree.codeml:
            runtime().debug("Importing codeml")
            self._codeml()
        else:
            runtime().debug("Already found codeml",self.family.alignments[0].tree.codeml.id)

        # Commit the session, close, and finish
        self.session.commit()
        self.session.close()

    def _index(self, name):
        n = name.split("#")[-1]
        if n.startswith("N"):
            n = n[1:]
        assert n.isdigit()
        return n

    def _tree(self):
        session = self.session

        # # Load the tree file and rename the taxa.
        # from Bio.Nexus.Nexus import Nexus
        # nex=Nexus(self.treeFile)
        # self.nexus = nex.trees[0]

        from Bio.Nexus.Trees import Tree as NewickTree
        tree_str = open(self.treeFile).read()
        self.nexus = NewickTree(tree_str)

        # Rename all the taxa.
        for id in self.nexus.get_terminals():
            node = self.nexus.node(id)
            node.data.taxon = self._index(node.data.taxon)

        # Create the DB object
        from hpf.hddb.db import Tree
        self.tree = Tree(alignment_key=self.alignment.id,
                         text=self.nexus.to_string(plain=False,plain_newick=True),
                         filename=self.treeFile
                         )
        session.add(self.tree)
        session.flush()

        # Now add in the node references
        self.nexus.name = self.tree.id
        assert self.tree.id != None
        runtime().debug("Added tree",self.tree)
        from hpf.hddb.db import TreeNodeFactory
        nodes = list(TreeNodeFactory().create(self.nexus))        
        for node in nodes:
            node.ancestor_node = node.ancestor.id if node.ancestor else None
            # This should add the new object into the session
            self.tree.nodes.append(node)
            #session.add(node)
            session.flush()
            
        runtime().debug("Appended",len(nodes),"tree nodes")
        session.flush()
        

        # Now import the diagnostic characters and reference the nodes.
        from hpf.amnh.oid import DiagCharsParser
        from hpf.hddb.db import TreeFactory
        biotree = TreeFactory(name_func=lambda node:str(node.id)).create(self.tree.nodes,self.tree.id)
        parser = DiagCharsParser(biotree)
        runtime().debug(self.treeDiagCharsFile)
        with open(self.treeDiagCharsFile) as handle:
            diagchars = list(parser.parse(handle))
            runtime().debug("DiagChars",len(diagchars))
            for d in diagchars:
                session.add(d)
        session.flush()

    def _codeml(self):
        if not self.codemlFile:
            return
        assert self.family.id != None
        assert self.tree.id != None

        # We need to convert the columns to the original alignment indices
        mapper = CulledColumnMapper(self.alignment, self.alignment.culled_columns)
        parser = PositiveSelectionParser()
        models = list(parser.parse(self.codemlFile))
        runtime().debug("Found",len(models),"models")
        for i,model in enumerate(models):
            model.tree_key = self.tree.id
            self.session.add(model)
            self.session.flush()
            ps = list(model.ps)
            runtime().debug("Found",len(ps),"sites in model",model.model)
            for j,site in enumerate(ps):
                site.codeml_key = model.id
                # Indices in CodeML start at 1, convert to 0 and then map
                orig = site.column
                site.column = mapper[site.column-1]
                runtime().debug("column",orig,"mapped to",site.column,site.probability)
                try:
                    self.session.add(site)
                except:
                    runtime().debug(i,":",j," failure on column",orig,"mapped to",site.column,site.probability)
                    raise
            runtime().debug("Finished with model")
            self.session.flush()
                
#        with open(self.codemlFile) as handle:
#            text = handle.read()
#        from hpf.hddb.db import CodeML
#        self.codeml = CodeML(tree_key=self.tree.id,
#                             filename=self.codemlFile,
#                             text=text)
#        self.session.add(self.codeml)
#        self.session.flush()
#        parser = LRTParser(self.alignment, self.alignment.culled_columns,self.codeml)
#        with open(self.codemlFile) as handle:
#            for selection in parser.parse(handle):
#                selection.codeml_key = self.codeml.id
#                self.session.merge(selection)
        runtime().debug("finished import codeml")
                
    def _alignment(self):
        session = self.session

        # Read the alignment
        from Bio import AlignIO
        with open(self.alignFile) as handle:
            align = AlignIO.read(handle,self.alignFormat)
        # Rename 'id' with the correct protein key
        for record in align:
            record.id = self._index(record.id)
        # Write to a text buffer and create the DB object
        text = StringIO()
        AlignIO.write([align],text,self.alignFormat)
        from hpf.hddb.db import Alignment
        self.alignment = Alignment(family_key=self.family.id,
                                   format=self.alignFormat,
                                   filename = self.alignFile,
                                   text=text.getvalue())
        # Add to session and flush
        session.add(self.alignment)
        session.flush()

        # Flip through the proteins in the alignment and add
        # the records.
        for record in align:
            protein_key = record.id
            assert protein_key!=0 and protein_key!=None, protein_key
            runtime().debug("protein: ",protein_key)
            from hpf.hddb.db import AlignmentProtein
            s = AlignmentProtein(alignment_key=self.alignment.id,
                                 protein_key=protein_key,
                                 sequence=str(record.seq))
            session.add(s)
            session.flush()
            
            # There may exist multiple alignments, but the definition
            # of membership in the family is done here.
            from hpf.hddb.db import FamilyProtein
            fs = FamilyProtein(family_key=self.family.id,
                               protein_key=protein_key,
                               seed=True)
            session.merge(fs)

        # Now read the colulmn culling log.  Indices start at 0 here.
        from hpf.hddb.db import AlignmentColcull, AlignmentSeqcull
        with open(self.alignColcullLog) as handle:
            for line in handle:
                column, gap, taxa, ratio = line.split()
                col = AlignmentColcull(alignment_key=self.alignment.id,
                                       column=column,
                                       gap_percentage=ratio)
                session.merge(col)
        with open(self.alignSeqcullLog) as handle:
            #rice#1182215    0.712765957446808
            for line in handle:
                parts= line.split()
                seq,score = parts
                seq = self._index(seq)
                #seq.split("#")[-1]
                if not seq.isdigit():
                    print parts,"SEQ:",seq
                    assert false
                cul = AlignmentSeqcull(alignment_key=self.alignment.id,
                                       protein_key=seq,
                                       score=score)
        session.flush()
        
    def _family(self):
        session = self.session
        from hpf.hddb.db import Family
        self.family = Family(name=self.familyName,
                             experiment_key=0)
        session.add(self.family)
        session.flush()

        
def main(*args):
    #pool = processor(synchronous=True)
    #runtime().debug("Using processor",pool)
    #pool.make_tasks(lambda: args)
    #consume(pool.run(import_dir))
    map(import_dir, args)

def import_dir(dir):
    assert os.path.isdir(dir)
    family_name = os.path.basename(os.path.abspath(dir))
    runtime().debug("Family Name",family_name)
    join = lambda *x: os.path.join(dir,*x)
    
        # Define the file names
    align_file = join("FAMILY.index")
    align_colcull_log = join("colcull.log")
    align_seqcull_log = join("seqcull.log")
    tree_file = join("oid.index.reroot")
    diagchars_file = join("diag.chars")

    from hpf.utilities import find
    codeml_file = join("%s_colcull.lrt" % family_name)
    print codeml_file, os.path.exists(codeml_file)
    codeml_file = codeml_file if os.path.exists(codeml_file) else None
    # Make sure all of the files exist
    for file in [tree_file,align_file,codeml_file,align_colcull_log,align_seqcull_log]:
        if file:
            assert os.path.exists(file),file 
    runtime().debug(family_name,
                    align_file,
                    tree_file,
                    codeml_file)
    oid = OIDImporter(familyName = family_name,
                      treeFile = tree_file,
                      treeDiagCharsFile = diagchars_file,
                      alignFile = align_file,
                      alignColcullLog = align_colcull_log,
                      alignSeqcullLog = align_seqcull_log,
                      codemlFile = codeml_file, 
                      alignFormat="fasta",
                      oid_key=oid_key)
    oid.merge()
        
#    pool = processor(synchronous=runtime().opt(SYNCHRONOUS))
#    runtime().debug("Using processor",pool)
#    pool.make_tasks(None)
#    consume(pool.run(None))
    
def _do(argv):
    r = runtime()
    r.description("""
    oid.py [-options] *oid-directory
    Import an OID family.
    """)
    #r.add_option(Flag(SYNCHRONOUS, "s", description="Run this script synchronously without any multi-processing"))
    args = r.parse_options(argv)
    main(*args)

if __name__=="__main__":
    _do(sys.argv[1:])
