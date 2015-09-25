import re
import sys, os
import fileinput
import datetime
import re

import os
FOLDER = os.path.abspath(os.path.dirname(__file__))
SCRIPTS_FOLDER = os.path.join(FOLDER,"scripts")
DATA_FOLDER = os.path.join(FOLDER,"data")
TOOLS_FOLDER = os.path.join(FOLDER,"tools")

from tempfile import NamedTemporaryFile, mkdtemp

# Coefficients are used as follows:
# resp = 
# c[0] + 
# c[1]*zscore +
# c[2]*prediction_contact_order +
# c[3]*convergence +
# c[4]*abs(log(length_ratio)) +
ALPHA_COEFFICIENTS = [-4.532031, 0.658800, 0.133027, -0.091581, -4.08231]
BETA_COEFFICIENTS = [-1.597025, 0.664228, 0.092968, -0.353935, -6.715978]
MIXED_COEFFICIENTS = [-4.095602, 0.673068, 0.051677, 0.025341, -5.160030]
# The probability is then 
# prob = 1/(1+1/exp(resp))


class McmNoTaskFactory(object):
    """A factory for creating McmData ORM objects (hpf.hddb.db.McmData) from MammothScore
    objects (hpf.mcm.mammoth.MammothScore). DOES NOT USE task chain structure (takes required
    SS and convergence values individually).
    """
    def __init__(self, mcmdb, percent_alpha, percent_beta, convergence_radius1, debug=False):
        """Parameters:
        mcmdb   - hpf.mcm.McmDB object containing list and data files with astral&scop information
        percent_alpha       - %alpha of the secondary structure prediction for decoy's sequence
        percent_beta        - %beta of the secondary structure prediction for decoy's sequence
        convergence_radius1 - radius1 value from convergence of clustering sequence's decoys
        ---
        NOTE that a created factory can correspond to only ONE (1) domain, for which the SS (and
        therefore %alpha and %beta), sequence length, and convergence values are all the same.
        This domain is the domain sent to Rosetta denovo structure pred, resulting in thousands
        of decoys that are then filtered and clustered to receive 25 cluster centers. The cluster
        centers are then run through mammoth, and all mammoth scores are MCM-scored with properties
        of the original single sequence.
        """
        self.debug = debug
        self.mcmdb = mcmdb
        self.percent_alpha = percent_alpha
        self.percent_beta  = percent_beta
        self.convergence_radius1 = convergence_radius1
        
    def create(self,score):
        """
        Create an McmData ORM object (hpf.hddb.db.McmData, hpf table mcm) from a Mammoth score.
        Entities:
          experiment_key  - structure ID of "experiment" source entry in McmDB obj that cluster center was compared against
          experiment      - McmDB Struct object with ID = experiment_key. All from mcmdb data.mammothDb file entries.
                          - Struct(id, sequence_key, length, percent_alpha/beta, sccs, astral_ac)
        Parameters:
          score           - hpf.mcm.mammoth.MammothScore object
        Returns an McmData ORM object corresponding to the given mammoth score
        """
        
        #?#TODO: Include values here to FULLY populate McmData ORM object. That way, return McmData 
        #?#TODO: objects can be added without modification to the DB
        ##NOTE: Leave outfile_mcm_result_key NULL (was used to reprsent MCM xml log files)
        # sequence_key
        # outfile_key
        # convergence_key
        # structure_key
        # astral_structure_key

        from hpf.hddb.db import McmData
        experiment_key = int(score.experiment.split(".")[0])
        experiment = self.mcmdb[experiment_key]

        prob,ratio,aratio,bratio,ssclass = self.probability(score=score, 
                                                            pred_percent_alpha=self.percent_alpha, 
                                                            pred_percent_beta=self.percent_beta, 
                                                            convergence=self.convergence_radius1, 
                                                            experiment=experiment)
        
        if self.debug: print "MCM probability: {0} {1} {2} {3} {4} {5}".format(experiment.id, 
                score.prediction, score.zscore, prob, ratio, ssclass)
        
        data= McmData(ratio = ratio,
                      aratio = aratio,
                      bratio = bratio,
                      percent_alpha = self.percent_alpha,
                      percent_beta  = self.percent_beta,
                      astral_percent_alpha = experiment.percent_alpha,
                      astral_percent_beta  = experiment.percent_beta,
                      sequence_length = int(score.num_p),
                      astral_sequence_length = int(score.num_e),
                      probability = prob,
                      sccs = experiment.sccs,
                      scop = self.mcmdb.scop,
                      )
        # "class" is a reserved word in python
        data.__dict__["class"] = ssclass
        # MammothScore object is saved in McmData object to make DB storage easier
        data.mammoth = score
        return data
            
    def probability(self, score, pred_percent_alpha, pred_percent_beta, convergence, experiment):
        """Calculates MCM stats on a given score with given prediction sequence properties
        Contains much magic from publications. Don't change the numbers...
        """    
        # Secondary structure class
        ssclass = 3
        if pred_percent_alpha >= 0.15 and pred_percent_beta < 0.15:
            ssclass = 1
        if pred_percent_alpha < 0.15 and pred_percent_beta >= 0.15:
            ssclass = 2
        
        # Length ratio
        ratio = float(score.num_p)/float(score.num_e) 
        
        # Alpha ratio
        if pred_percent_alpha != 0 and float(experiment.percent_alpha) != 0:
            aratio = float(pred_percent_alpha) / float(experiment.percent_alpha)
        elif pred_percent_alpha == 0 and float(experiment.percent_alpha) == 0:
            aratio = 1.0
        elif pred_percent_alpha == 0:
            aratio = 0.0067379
        elif float(experiment.percent_alpha) == 0:
            aratio = 148.41
        else:
            assert False, "Problem with percent_alpha. This cannot happen..."
        
        # Beta ratio
        if pred_percent_beta != 0 and float(experiment.percent_beta) != 0:
            bratio = float(pred_percent_beta) / float(experiment.percent_beta)
        elif pred_percent_beta == 0 and float(experiment.percent_beta) == 0:
            bratio = 1.0
        elif pred_percent_beta == 0:
            bratio = 0.0067379
        elif float(experiment.percent_beta) == 0:
            bratio = 148.41
        else:
            assert False, "Problem with percent_beta. This cannot happen..."

        # Coefficients
        from numpy import log,abs,sum,exp
        variables = [1.0, 
                     score.zscore, 
                     score.prediction_contact_order, 
                     convergence, 
                     abs(log(ratio))]
        
        if ssclass==1:
            coeff = ALPHA_COEFFICIENTS
        elif ssclass==2:
            coeff = BETA_COEFFICIENTS
        elif ssclass==3:
            coeff = MIXED_COEFFICIENTS
        else:
            assert False, "Problem with SS class. This cannot happen..."

        assert len(variables)==len(coeff), "length of vars and coefficients should be equal"
        weighted = [float(v)*co for v,co in zip(variables,coeff)]
        resp = sum(weighted)
        prob = 1.0/(1.0 + 1.0/exp(resp))
        assert prob>=0 and prob <=1.0, "This is supposed to be a probability! (between 0 and 1) %f" % prob
        return prob,ratio,aratio,bratio,ssclass
        

class McmFactory(object):
    
    def __init__(self, convergence, mcmdb):
        self.convergence = convergence
        self.mcmdb = mcmdb
        
    def create(self,score):
        """
        Create an McmData ORM object (hpf.hddb.db.McmData, hpf table mcm) from a Mammoth score.
        Patrick:
          Conspicuously missing from here is an 'outfile_key'. This is tricky.  McmData objects 
          store EVERYTHING compiled from Rosetta clustering runs, mammoth, and the mcm database.
        Entities:
          experiment_key  - structure ID of "experiment" source entry in McmDB obj that cluster center was compared against
          experiment      - McmDB Struct object with ID = experiment_key. All from mcmdb data.mammothDb file entries.
                          - Struct(id, sequence_key, length, percent_alpha/beta, sccs, astral_ac)
          prediction      - hpf.mcm.cluster.ClusterCenter object from convergence task (TODO: CHANGE)
        Parameters:
          score           - hpf.mcm.mammoth.MammothScore object
        Returns an McmData ORM object corresponding to the given mammoth score
        """
        from hpf.hddb.db import McmData
        experiment_key = int(score.experiment.split(".")[0])
        experiment = self.mcmdb[experiment_key]
        prediction = self.convergence[score.prediction]

        prob,ratio,aratio,bratio,ssclass = self.probability(score,prediction,experiment)
        e = experiment
        p = prediction
        s = score
        print "probability: ", e.id, p.index, s.zscore, prob, ratio, ssclass
        
        data= McmData(ratio = ratio,
                      aratio = aratio,
                      bratio = bratio,
                      percent_alpha = p.percent_alpha,
                      percent_beta = p.percent_beta,
                      astral_percent_alpha = e.percent_alpha,
                      astral_percent_beta = e.percent_beta,
                      sequence_length = s.num_p,
                      astral_sequence_length = int(s.num_e),
                      probability = prob,
                      sccs = e.sccs,
                      scop = self.mcmdb.scop,
                      )
        # "class" is a reserved word in python
        data.__dict__["class"] = ssclass
        data.mammoth = score
        return data
            
    def probability(self, score, prediction, experiment):
        pred = prediction
        
        # Secondary structure class
        ssclass = 3
        if pred.percent_alpha >= 0.15 and pred.percent_beta < 0.15:
            ssclass = 1
        if pred.percent_alpha < 0.15 and pred.percent_beta >= 0.15:
            ssclass = 2
        
        # Length ratio
        ratio = float(score.num_p)/float(score.num_e) 
        
        # Alpha ratio
        if pred.percent_alpha != 0 and float(experiment.percent_alpha) != 0:
            aratio = float(pred.percent_alpha) / float(experiment.percent_alpha)
        elif pred.percent_alpha == 0 and float(experiment.percent_alpha) == 0:
            aratio = 1.0
        elif pred.percent_alpha == 0:
            aratio = 0.0067379
        elif float(experiment.percent_alpha) == 0:
            aratio = 148.41
        else:
            assert False, "Problem with percent_alpha. This cannot happen..."
        
        # Beta ratio
        if pred.percent_beta != 0 and float(experiment.percent_beta) != 0:
            bratio = float(pred.percent_beta) / float(experiment.percent_beta)
        elif pred.percent_beta == 0 and float(experiment.percent_beta) == 0:
            bratio = 1.0
        elif pred.percent_beta == 0:
            bratio = 0.0067379
        elif float(experiment.percent_beta) == 0:
            bratio = 148.41
        else:
            assert False, "Problem with percent_beta. This cannot happen..."

        # Coefficients
        from numpy import log,abs,sum,exp
        variables = [1.0, 
                     score.zscore, 
                     score.prediction_contact_order, 
                     pred.convergence, 
                     abs(log(ratio))]
        
        if ssclass==1:
            coeff = ALPHA_COEFFICIENTS
        elif ssclass==2:
            coeff = BETA_COEFFICIENTS
        elif ssclass==3:
            coeff = MIXED_COEFFICIENTS
        else:
            assert False, "this cannot happen..."

        assert len(variables)==len(coeff), "length of vars and coefficients should be equal"
        weighted = [float(v)*co for v,co in zip(variables,coeff)]
        resp = sum(weighted)
        prob = 1.0/(1.0 + 1.0/exp(resp))
        assert prob>=0 and prob <=1.0, "This is supposed to be a probability! %f" % prob
        return prob,ratio,aratio,bratio,ssclass
        
        
class McmDB(object):
    """
    Load astral id's or database id's into a set. For entry fields, see the .add() method.
    """
    def __init__(self,list_file=None,data_file=None,scop=1.75):
        self._ids = set()
        self._dict = dict()
        self.list_file = list_file
        self.data_file = data_file
        self.scop = scop
    
    @staticmethod
    def load(data_file,db=None):
        """
        Static method to load a new database from a data file.
        @param db: Alternatively add entries to existing db.
        """
        file = data_file
        db = McmDB()
        # Pull the header line from 'file'
        keys = file.readline().strip().split()
        class Struct(object):
            pass
        # For remaining data lines in 'file', create Struct object with attributes from file header set to the line's values
        # EG: Struct(structure_key, sequence_key, length, percent_alpha, percent_beta, sccs, astral_ac) and set Struct.id=structure_key
        # Then add the entry to the McmDB object (id to self._ids set, entry to self._dict {id -> Struct})
        for line in file:
            entry = Struct()
            for value, key in zip(line.strip().split(),keys):
                setattr(entry, key, value)
            entry.id = entry.structure_key
            db.add_entry(entry)
        return db
            
    def add_entry(self, entry):
        self._ids.add(entry.id)
        self._dict[entry.id]=entry
    
    def add(self, 
            id, 
            sequence_key=None, 
            length=None,
            percent_alpha=None,
            percent_beta=None,
            sccs=None,
            astral_ac=None):
        """
        Add an id to this databases definition set
        'id' is structure key in data.mammothDb file
        """
        class Struct(object):
            pass
        entry = Struct()
        entry.id = id
        entry.sequence_key = sequence_key
        entry.length = length
        entry.percent_alpha = percent_alpha
        entry.percent_beta = percent_beta
        entry.sccs = sccs
        entry.astral_ac = astral_ac
        return self.add_entry(entry)

    def __getitem__(self, key):
        return self._dict[str(key)]

    def add_astral(self, astral_id, session):
        """
        Add an astral id, searching the database to make sure it exists.
        @return: tuple: (self.add(id), id) : id corresponding to the db astral record
        """
        from hpf.hddb.db import Astral
        astral = session.query(Astral).filter(Astral.stype+
                                                    Astral.pdbid+
                                                    Astral.part == 
                                                    astral_id).first()
        assert astral!=None
        id = astral.id
        return (self.add(id),id)
        
class McmDBExporter(object):
    """
    Writes an MCM data file, calculating SS values using DSSP.
    """
    
    # Order of columns for data file
    columns = ["structure_key",
               "sequence_key",
               "length", 
               "percent_alpha", 
               "percent_beta", 
               "sccs", 
               "astral_ac"]
    
    def __init__(self, 
                 mcmdb, 
                 mammoth_list="list.mammoth", 
                 info_file="data.mammothDb", 
                 dir=None):
        self.mcmdb = mcmdb
        self.mammoth_list = mammoth_list
        self.info_file = info_file
        self.dir = dir
        
    def __enter__(self):
        if not self.dir:
            self.dir = mkdtemp()
        from hpf.hddb.db import Session
        self.session = Session()
        self.list_handle = open(self.mammoth_list,"w")
        self.info_handle = open(self.info_file,"w")
        
        print >>self.info_handle, "\t".join(McmDBExporter.columns)
        print >>self.list_handle, "MAMMOTH List\n%s" % self.dir
        return self
        
    def __exit__(self, type, value, traceback):
        self.session.close()
        self.session = None
        self.list_handle.close()
        self.info_handle.close()

    def export(self):
        for id in self.mcmdb._ids:
            try:
                self._write(self._dict(id))
            except:
                print "error on",id
                continue

    def _write(self, d):
        data = [d[key] for key in McmDBExporter.columns]
        astral = d["astral"]
        print >>self.info_handle, "\t".join([str(value) for value in data])
        f = "%i.pdb" % astral.structure_key
        print >>self.list_handle, f
        full = os.path.join(self.dir,f)
        with open(full,"w") as handle:
            handle.write(astral.structure.text)
        
        
    def _dict(self, id):
        """
        @return:  
        """
        values = {}
        from hpf.hddb.db import Astral
        astral = self.session.query(Astral).get(id)
        dssp = self._dssp(astral.structure)
        
        ss = []
        first_model = list(astral.structure.pdb)[0]
        for i,res in enumerate(first_model.get_residues()):
            chain = res.get_parent().get_id() 
            hetero,seq,insertion = res.get_id()
            key = (chain,(' ',seq,insertion))
            try:
                aa, s, accessibility = dssp[key]
                ss.append(s)
            except:
                continue
        length = i+1
        assert float(abs(length-len(ss)))/float(length)<0.1
        alpha = [a for a in ss if a in ('H','G','I')]
        beta = [a for a in ss if a in ('E','B')]
        
        values['percent_alpha'] = float(len(alpha))/float(length)
        values['percent_beta'] = float(len(beta))/float(length)
        values["astral"] = astral
        values["length"] = length
        values["sccs"] = astral.sccs
        values["astral_ac"] = astral.stype+astral.pdbid+astral.part
        values["sequence_key"] = astral.sequence_key
        values["structure_key"] = astral.structure_key
        
        return values
            
    def _dssp(self, structure):
        structure_temp = NamedTemporaryFile(dir = self.dir)
        structure_file = structure_temp.name
        with open(structure_file,"w") as handle:
            handle.write(structure.text)
        
        dssp_temp = NamedTemporaryFile(dir = self.dir)
        cmd = "dssp %s > %s" % (structure_file, dssp_temp.name)
        import subprocess
        print cmd
        subprocess.check_call(cmd,shell=True)
        from Bio.PDB.DSSP import make_dssp_dict
        dict, keys = make_dssp_dict(dssp_temp.name)
        # This can be dangerous...
        #os.system("rm "+out_file)
        return dict
        
