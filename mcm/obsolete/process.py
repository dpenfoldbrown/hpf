'''
Created on Apr 30, 2010

@author: Patrick
'''
import sys
import os
from hpf.utilities.process import TaskHandler,OutputTask, Task        
from hpf.runtime import runtime
from hpf.pdb.mammoth import Mammoth
from hpf.mcm import McmFactory
from hpf.mcm.mammoth import MastodonCL,MammothList,MammothMultiParser
from hpf.mcm.cluster import ClusterCenter

class McmReRun(TaskHandler):
    """
    Re-run MCM for a given set of sequence keys.
    """
    
    def __init__(self, dir, db, sequence_key,**kwargs):
        super(McmReRun,self).__init__(**kwargs)
        self.db = db
        self.sequence_key = sequence_key
        self.convergence_task = XMLConvergenceTask(self.sequence_key)
        self.add_task(self.convergence_task)
        self.mammoth_task = MammothTask(self.convergence_task,
                                        self.db)
        self.add_task(self.mammoth_task)
        self.mcm_task = McmTask(self.convergence_task,
                                self.mammoth_task,
                                self.db)
        self.add_task(self.mcm_task)
        self.db_task = McmDbTask(self.sequence_key,
                                 self.db,
                                 self.convergence_task,
                                 self.mammoth_task,
                                 self.mcm_task)
        self.add_task(self.db_task)

class McmDbTask(Task):
    """
    Handle the saving of MCM and convergence records.
    """
    
    def __init__(self, sequence_key, mcmdb, convergence, mammoth, mcm, *args, **kwargs):
        super(McmDbTask,self).__init__(*args,**kwargs)
        self.sequence_key = sequence_key
        self.mcmdb = mcmdb
        self.convergence = convergence
        self.mammoth = mammoth
        self.mcm = mcm
        
    def _convergence(self):
        """
        Upload the convergence values, cluster centers, and decoy structures
            as necessary.
        @note: This is using some complicated sqlalchemy relations and care
            should be taken to make sure the relationships are being maintained.
        """
        from hpf.hddb.db import RosettaConvergence, RosettaCluster
        from sqlalchemy.sql import and_
        outfile_key = self.convergence.outfile_key
        rb = self.convergence.robetta
        cv = self.session.query(RosettaConvergence).filter(RosettaConvergence.outfile_key==outfile_key).first()
        cv = cv if cv != None else RosettaConvergence()
        
        # Update the convergence record
        from hpf.hddb.db import BaseStruct
        vars = {
            'outfile_key': outfile_key,
            'target':os.path.basename(self.convergence.params.outfile),
            'radius1':rb.radius1,
            'size1':rb.size1,
            'radius2':rb.radius2,
            'size2':rb.size2,
            'total_decoys':rb.total_n_decoys
            }
        cv._update_vars(vars)
        assert cv.outfile_key != None, cv
        self.session.merge(cv)
        self.cv = cv
        self.session.flush()
        runtime().debug("Convergence record",cv)
        assert cv.id!=None, cv
        
        # index the saved cluster centers for getting mammoth scores
        self.structures = {}
        for cluster_center in set(self.convergence.structures.values()):
            cc = self.session.query(RosettaCluster).filter(and_(RosettaCluster.convergence_key==cv.id,
                                                                RosettaCluster.index==cluster_center.index
                                                                )).first()
            cc = cc if cc != None else RosettaCluster()
            vars = {
                'index': cluster_center.index,
                'size': cluster_center.size,
                'rank': cluster_center.rank,
                'convergence_key':self.cv.id
                }
            cc._update_vars(vars)
            # Upload the decoy structure if it doesn't already exist
            if cc.structure == None or cc.structure_key==None:
                from hpf.hddb.db import Structure
                structure = Structure(sequence_key=self.sequence_key,
                                      structure_type="decoy",
                                      comment=self.convergence[cc.index].decoy_file,
                                      text=self.convergence.centers[cc.index].decoy()
                                      )
                assert structure.sequence_key!=None,structure
                self.session.add(structure)
                self.session.flush()
                cc.structure_key = structure.id
                cc.structure = structure
                runtime().debug("\tnew Structure",structure)
            else:
                runtime().debug("\texisting Structure",cc.structure_key,cc.structure)


            assert cc.convergence_key!=None and cc.structure_key!=None,cc
            self.session.merge(cc)
            self.session.flush()
            runtime().debug("\tCluster record",cc)
            self.structures[cc.index] = cc.structure
            
    def _mcm(self):
        """Save the top 5 MCM records."""
        from hpf.hddb.db import MammothFactory
        factory = MammothFactory()
        # Returns the 5 objects with greatest 'probability' values. Objects are hpf.hddb.db.McmData objs
        best = sorted(self.mcm,reverse=True)[0:5]
        for mcmdata in best:
            cluster_index = int(mcmdata.mammoth.prediction.split(".")[0][6:])
            structure_key = self.structures[cluster_index].id
            astral_structure_key = int(mcmdata.mammoth.experiment.split(".")[0])
            
            score = factory.create(mcmdata.mammoth)
            score.p_structure_key = structure_key
            score.e_structure_key = astral_structure_key
            
            mcmdata.sequence_key = self.sequence_key
            mcmdata.outfile_key = self.cv.outfile_key
            mcmdata.convergence_key = self.cv.id
            mcmdata.structure_key = structure_key
            mcmdata.astral_structure_key = astral_structure_key
            
            self.session.merge(score)
            self.session.merge(mcmdata)
            self.session.flush()
            #assert score.p_structure_key!=None and score.e_structure_key!=None, score
            
        #decoy_9057.pdb  987697.pdb
        
    def _pre(self):
    # Sets up DB session
        from hpf.hddb.db import Session, tunnel, rebind, url
        from sqlalchemy.orm import sessionmaker
        from sqlalchemy import create_engine
        tunnel()
        # global application scope.  create Session class, engine
        Session = sessionmaker()
        engine = create_engine(url+'hpf',pool_recycle=300,echo=True,pool_timeout=10)
        # local scope, such as within a controller function
        # connect to the database
        connection = engine.connect()
        # bind an individual Session to the connection
        self.session = Session(bind=connection)
        
    def _post(self):
    # Flushes and closes session
        self.session.commit()
        self.session.close()
    
    def _do(self):
        self._convergence()
        self._mcm()
        
        
class McmTask(Task):
    """
    Prints the list of scores that is the MammothTask obj.
    Runs MammothTask scores through McmFactory's create() method to produce McmData objs (1 per score)
    Treated as an iterable, iterates over mcmdata, the list of McmData objects from scores.
    """

    def __init__(self, convergence_task, mammoth_task, mcmdb, **kwargs):
        super(McmTask,self).__init__(**kwargs)
        self.convergence_task = convergence_task
        self.mammoth_task = mammoth_task
        self.mcmdb = mcmdb

    def _do(self):
        from itertools import imap
        print "McmTask:: mammoth_task.scores (pred, exp, zscore):"
        for score in self.mammoth_task.scores:
            print score.prediction, score.experiment, score.zscore
        
        factory = McmFactory(self.convergence_task, self.mcmdb)
        self.mcmdata = imap(factory.create,self.mammoth_task.scores)
        
    def __iter__(self):
        return iter(self.mcmdata)
        
class MammothTask(OutputTask):
    """
    Creates a mammoth list from the set of decoys and runs mammoth against
        the given MCM database.
    Treated as an iterable, iterates over self.scores
    """
    
    def __init__(self, decoys, mcmdb):
        super(MammothTask,self).__init__("data.mammoth")
        self.decoys = decoys
        self.mcmdb = mcmdb
    
    def _passed(self):
        print "_passed"
        with open("data.mammoth") as handle:
            for line in handle:pass
            if -1 == line.find("NORMAL_EXIT"):
                raise Exception, "Mammoth did not exit normally (NORMAL_EXIT not found)"
            else:
                print "Mammoth exited normally"
        with open("data.mammoth") as handle:
            self.scores = list(MammothMultiParser().parse(handle))

    def _check(self):
    with open("data.mammoth") as handle:
        import re
        #kdrew: tests to see if score lines in data.mammoth (ie. lines starting with ':') match scores in data structure
        if len(list(line for line in handle if re.match("^:",line))) != len(self.scores):
            raise Exception, "Not all scores from data.mammoth were read"
        else:
            print len(self.scores), " scores were read"

    def _post(self):
        print "_post"
        runtime().debug("Mammoth Scores:",len(self.scores))
        self._check()

    def _do(self):
        print "_do"
        # MammothList(file, directory, *structures)
        self.mammoth_list = MammothList("cmd.mammoth", os.getcwd(), *self.decoys)
        
        from hpf.hddb.mcm import TOOLS_FOLDER
        command = os.path.join(TOOLS_FOLDER,"mammoth")
        
        with self.mammoth_list as list_file:
            runtime().debug("Mammoth list file",list_file)
            # Create a Mastodon (Mammoth + contact order) command line string
            self.cl = MastodonCL(contact_order=True,
                            experiment=self.mcmdb.list_file,
                            prediction=list_file,
                            cwd=os.getcwd(),
                            output="data.mammoth",
                            level=0,
                            verbose=0,
                            command = command)
            # Run Mammoth with the custom (w/contact order) cmdline string
            Mammoth(self.cl, parse=False).run()
        
        # Open the output file. Check for successful exit, then parse scores
        with open(self.cl.output) as handle:
            for line in handle:pass
            if -1 == line.find("NORMAL_EXIT"):
                raise Exception, "Mammoth did not NORMAL_EXIT"
            else:
                print "NORMAL_EXIT found"
        with open(self.cl.output) as handle:
            self.scores = list(MammothMultiParser().parse(handle))
        print self.scores

    def __iter__(self):
        return iter(self.scores)


class XMLConvergenceTask(OutputTask):
    """
    Gather convergence values and structures from an MCM XML result file for
        a set of sequence keys.  Supports iteration over decoy filenames
        and dict like access for ClusterCenter convergence info by filename
        or structure_key.
    """

    def __init__(self, sequence_key):
        super(XMLConvergenceTask,self).__init__()
        self.sequence_key = sequence_key
        self.structures = {}
        self.centers = {}
        self.session = None
        self.log_file = None
    
    def _xml_log(self):
    """ Gets MCM clustering XML outfile from database by self sequence_key
        Parses outfile and returns McmLogFile object (mcm/xml.py)
    """
        if self.log_file:
            runtime().debug("Found logfile")
            return self.log_file
        from hpf.mcm.xml import McmLogFile
        from hpf.hddb.db import McmResultFile, Session
        self.session = Session()
        # For each sequence grab the most recent XML file
        sequence_key = self.sequence_key
        runtime().debug("Querying")
        result_file = self.session.query(McmResultFile
                           ).filter(McmResultFile.sequence_key==sequence_key
                                    ).order_by(McmResultFile.scop.desc()
                                               ).first()
        if not result_file:
            raise Exception("No MCM result file for sequence {0} found".format(sequence_key))
        runtime().debug("Found mcm result file for sequence {0}".format(sequence_key),result_file)
        self.outfile_key = result_file.outfile_key
        log_file = McmLogFile.parseString(result_file.text)
        return log_file

    def _pre(self):
        """Gather all cluster centers from the latest XML file for self sequence_key"""
        from itertools import islice
        
        log_file = self._xml_log()
        params = log_file.params()
        self.params = params
        
        convergence = log_file.convergence()
        self.robetta = convergence
        
        length = len(params.sequence)
        percent_alpha = float(params.ss.count('H'))/float(length)
        percent_beta = float(params.ss.count('E'))/float(length)
       
        for center in log_file.cluster_centers():
            # McmLogFile.cluster_centers is dict {cluster_index => xml.ClusterCenter obj}
            runtime().debug("Center",center.index)
            
            #TODO: Change to get seq length, ss content, etc from first mcmdata entry for sequence
            #entry = list(islice(center.entries(),1))[0]
            c = ClusterCenter("decoy_%i.pdb" % center.index,    # decoy_file
                              float(convergence.radius1),       # convergence
                              center.size,                      # size
                              center.index,                     # index
                              center.rank,                      # rank
                              int(convergence.total_n_decoys),  # n_decoys_in_outfile
                              os.path.basename(params.outfile), # target (aka silent file)
                              self.sequence_key,                # sequence_key
                              length,                           # sequence_length
                              percent_alpha,                    # percent_alpha
                              percent_beta                      # percent_beta
                )
            # Add to structures dict by decoy_filename and index
            self.structures[c.decoy_file] = c
            self.structures[c.index] = c
            
            # Index the XML entries so we can grab decoy atom records
            self.centers[c.index] = center
            runtime().debug("\tDecoy for cluster center created",c.decoy_file)
        
        # Files list: all unique decoy filenames in structures dict
        files = list(set([c.decoy_file for c in self.structures.values()]))
        runtime().debug("Need",files)
        
        # Add files list to required existing output for task
        self.append_output(*files)
        self.decoy_files = files

    def _export(self,id):
        """Write a structure with the given id to the 'decoy_file'"""
        from hpf.hddb.db import Structure
        filename = self[id].decoy_file
        # xml.ClusterCenter.decoy() gets atomrecord from XML object for center
        decoy = self.centers[id].decoy()
        with open(filename,"w") as handle:
            print >>handle, decoy
        return filename

    def _post(self):
        if self.session:
            self.session.close()
    
    def _do(self):
        """Export all parsed cluster centers to decoy files"""
        cluster_indices = self.centers.keys() 
        self.decoy_files = map(self._export,cluster_indices)
    
    def __iter__(self):
        """Iterator over the decoy files"""
        return iter(self.decoy_files)
    
    def __getitem__(self, key):
        """Dict-like access for decoy names and structure keys."""
        return self.structures[key]
    
class ConvergenceTask(OutputTask):
    """
    Gather convergence values and structures from McmData for
        a set of sequence keys.  Supports iteration over decoy filenames
        and dict like access for ClusterCenter convergence info by filename
        or structure_key.
    """

    def __init__(self, *sequence_keys):
        assert False, "Don't use this."
        super(ConvergenceTask,self).__init__()
        self.sequence_keys = sequence_keys
        self.structures = {}
    
    def _pre(self):
        """
        Gather all structures from the last scop run on these sequences. 
        """
        from hpf.hddb.db import McmData, Session
        from sqlalchemy.sql import func, and_
        self.session = Session()
        # Get all McmData objs with latest SCOP from given seq keys, make ClusterCenters from them, and add to structures dict (by key and filename)
        # Mad redundant querying
        
        ## NOTE: This will not work as it stands. The McmData DB ORM object has been re-mapped to the mcm table
        # (instead of the mcmData table, which has the fields required for ClusterCenter.create()'s parameter but is not written)
        
        results = self.session.query(McmData.sequence_key, func.max(McmData.scop)).group_by(McmData.sequence_key, \
                McmData.scop).filter(McmData.sequence_key.in_(self.sequence_keys)).distinct().all()
        for sequence_key, scop in results:
            mcmdata = self.session.query(McmData).filter(and_(McmData.sequence_key == sequence_key, \
                    McmData.scop == scop)).all()
            for data in mcmdata: 
                c = ClusterCenter.create(data)
                self.structures[data.structure_key] = c
                self.structures[data.prediction_file] = c   

        #runtime().debug("Found",self.structures.values)
        files = list(set([c.decoy_file for c in self.structures.values()]))
        runtime().debug("Need",files)
        self.append_output(*files)
        self.decoy_files = files

    def _export(self,id):
        """Write a structure with the given id to the 'decoy_file'"""
        from hpf.hddb.db import Structure
        filename = self[id].decoy_file
        structure = self.session.query(Structure).get(id)
        with open(filename,"w") as handle:
            print >>handle, structure.text
        return filename

    def _post(self):
        self.session.close()
    
    def _do(self):
    ''' Export structure files from all structure keys
    '''
        structure_keys = list(set([c.structure_key for c in self.structures.values()]))
        self.decoy_files = map(self._export,structure_keys)
    
    def __iter__(self):
        """Iterator over the decoy files"""
        return iter(self.decoy_files)
    
    def __getitem__(self, key):
        """Dict like access for decoy names and structure keys."""
        return self.structures[key]
    
class ClusterTask(OutputTask):
    
    def __init__(self, outfile):
        pass
