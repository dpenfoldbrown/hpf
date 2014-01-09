'''
Necessary clustering tools for MCM.
Created on Apr 28, 2010
@author: Patrick
'''

## Example robetta_clusterer cmd file:
#OUTPUT_FILE data.cluster
#TARGET /scratch/pcw216/tmp/mcm/sandbox/jobs/au339/au339.15000.10000.out PAGSEGYVDGFDHTAWRYLMSPYISAYKLGLSEPYINFESLFYWYRPTPKSATATADSL


## Example header of robetta_clusterer outfile:
#THRESHOLD: 3.683921 TOP_CLUSTER_SIZE: 279
#standard_thresholds: size1= 186 threshold1= 3.342081 size2= 465 threshold2= 4.056837 total_decoys= 9306

import os
import re
from hpf.utilities.temp import ControlFile
from hpf.utilities.process import Process
#
#CREATE TABLE `mcmData` (
#  `id` int(11) NOT NULL AUTO_INCREMENT,
#  `sequence_key` int(11) NOT NULL DEFAULT '0',
#  `outfile_key` int(11) NOT NULL DEFAULT '0',
#  `structure_key` int(11) NOT NULL,
#  `n_decoys_in_outfile` int(11) NOT NULL DEFAULT '0',
#  `cluster_center_rank` int(11) NOT NULL DEFAULT '0',
#  `target` varchar(50) NOT NULL DEFAULT '',
#  `convergence` double NOT NULL DEFAULT '0',
#  `prediction_file` varchar(150) NOT NULL DEFAULT '',
#  `decoy_name` varchar(20) NOT NULL DEFAULT '',
#  `prediction_percent_alpha` double NOT NULL DEFAULT '0',
#  `prediction_percent_beta` double NOT NULL DEFAULT '0',
#  `prediction_sequence_length` int(11) NOT NULL DEFAULT '0',
#  `cluster_center_size` int(11) NOT NULL DEFAULT '0',
#  `cluster_center_index` int(11) NOT NULL DEFAULT '0',
#  `timestamp` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
#  PRIMARY KEY (`id`),
#  UNIQUE KEY `outfile_key` (`outfile_key`,`decoy_name`,`experiment_file`),
#  KEY `sequence_key` (`sequence_key`),
#  KEY `outfile_key` (`outfile_key`),
#  KEY `structure_key` (`structure_key`),

class ClusterCenter(object):
    """CLASS DEPRECATED"""

    def __init__(self,
                 decoy_file,
                 convergence,
                 size,
                 index,
                 rank,
                 n_decoys_in_outfile,
                 target,
                 sequence_key,
                 sequence_length,
                 percent_alpha,
                 percent_beta,
                 outfile_key=None,
                 structure_key=None
                 ):
        
        raise Warning("hpf.mcm.cluster.ClusterCenter class is currently DEPRECATED")
        
        self.decoy_file = decoy_file
        self.convergence = convergence
        self.size = size
        self.index = index
        self.rank = rank
        self.n_decoys_in_outfile = n_decoys_in_outfile
        self.target = target
        self.sequence_key = sequence_key
        self.sequence_length = sequence_length
        self.percent_alpha = percent_alpha
        self.percent_beta = percent_beta
        self.convergence = convergence
        self.outfile_key=outfile_key
        self.structure_key = structure_key
    
    @staticmethod
    def create(mcmdata):
        """
        Takes an McmData DB ORM object (representing hpf.mcm table)
        
        NOTE: This will not currently work. McmData has been mapped to hpf.mcm in db.py ORM,
        where this requires the fields from the hpf.mcmData table (for which no ORM object
        is written.
        """
        data = mcmdata
        return ClusterCenter(data.prediction_file, 
                             data.convergence, 
                             data.cluster_center_size, 
                             data.cluster_center_index, 
                             data.cluster_center_rank, 
                             data.n_decoys_in_outfile,
                             data.target,
                             data.sequence_key,
                             data.prediction_sequence_length,
                             data.prediction_percent_alpha,
                             data.prediction_percent_beta,
                             data.outfile_key,
                             data.structure_key)


class RobettaClusterCenter():
    """A class to store a Robetta cluster center. List of these provided to RobettaCluster obj.
    Params:
        rank    - the rank (0->best) of the cluster
        index   - the clusterer-internal index of the cluster center decoy
        rosetta_id  - the S_XXXX_XXXX rosetta "description" given in denovo result file
        size    - the size of the cluster
        pdb_file- OPTIONAL filename of pdb file that holds atom record of cluster center
    """
    def __init__(self, rank, index, rosetta_id, size, pdb_file=None):
        self.rank = rank
        self.index = index
        self.rosetta_id = rosetta_id
        self.size = size
        self.pdb_file = pdb_file

    def __repr__(self, ):
        return "<RobettaClusterCenter: rank {0}, index {1}, id {2}, size {3}>".format(self.rank, self.index, self.rosetta_id, self.size)

    def get_atom_record(self, ):
    # Returns the atom record (whole pdb file as a string) of the center.
    # NOTE: If pdb_file is not defined, returns None
        if not self.pdb_file:
            return None
        if not os.path.isfile(self.pdb_file):
            raise Exception("{0} pdb file '{1}' is not a valid file".format(self, self.pdb_file))
        with open(self.pdb_file) as handle:
            atom_record = handle.read()
        return atom_record

class RobettaConvergence():
    """A class to store convergence values parsed from Robetta clusterer result file. Passed to RobettaCluster obj.
    Params:
        result_file - the clusterer result file parsed to get these values
        size1   - the "size1" value from clusterer result file
        radius1 - ..  "threshold1" ..
        size2   - ..  "size2"      ..
        radius2 - ..  "threshold2" ..
        total_decoys    - the total number of decoys clustered, parsed from result file
    """
    def __init__(self, result_file, size1, radius1, total_decoys, size2=None, radius2=None):
        self.result_file  = result_file
        self.size1        = size1
        self.radius1      = radius1
        self.size2        = size2
        self.radius2      = radius2
        self.total_decoys = total_decoys

    def __repr__(self, ):
        return "<RobettaConvergence: from {0}, size1 {1}, radius1 {2}, total decoys {3}>".format(self.result_file, self.size1, self.radius1, self.total_decoys)

class RobettaCluster():
    """A class to store information parsed from Robetta clusterer outfile
    Params:
        result_file - the filename of the clusterer result file parsed for these values
        threshold   - the inital threshold ('THRESHOLD') value at the top of the cluster file
        convergence - a RobettaConvergence object, storing radiuss and sizes (1 per Clusterer run/RobettaCluster object)
        centers     - a list of RobettaClusterCenter objects representing cluster centers from the result file
    Treated as an iterable, iterates over self.centers list. 
    Treated as a hashable, returns center obj from self.centers_dict (key: center index OR center rosetta ID)
    """
    def __init__(self, result_file, threshold, top_cluster_size, convergence=None, centers=None):
        self.result_file      = result_file
        self.threshold        = threshold
        self.top_cluster_size = top_cluster_size
        self.convergence      = convergence
        self.centers          = centers
        self.centers_dict     = None

        if self.centers:
            self.centers_dict = self._create_centers_dict()

    def _create_centers_dict(self, ):
    # Create dict (int) index/'rosetta_id' => RobettaClusterCenter obj
        centers_dict = dict()
        for center in self.centers:
            centers_dict[center.index] = center
            centers_dict[center.rosetta_id] = center
        return centers_dict

    def __repr__(self, ):
        repr = "<RobettaCluster: from {0}".format(self.result_file)
        if self.centers:
            repr += ", centers: {0}".format(len(self.centers))
        if self.convergence:
            repr += ", convergence: {0}".format(self.convergence.radius1)
        repr += ">"
        return repr

    def __iter__(self, ):
    # Treated as an iterable, return iter over centers list (if available)
        if self.centers:
            return iter(self.centers)
        else:
            return None

    def __getitem__(self, key):
    # Treated as a hashable, return value from centers_dict (key is center.index or center.rosetta_id)
        if self.centers_dict:
            return self.centers_dict[key]
        else:
            return None


class RobettaClusterer():
    """A class to wrap the whole functionality of running the robetta_cluster program
    Parameters:
        command_file- filename to populate with clusterer parameters, described below
        decoy_file  - the file (SILENT format) of decoys to be clustered
        sequence    - the sequence of the decoys to be clustered (if None, parsed from decoy_file)
        outfile     - the file clusterer will output results to (added to 'control_file')
        log_file    - optional file for keeping clusterer STDOUT
        executable  - the clusterer binary
        path        - an optional path. If given, will create file and run there (and read flat filenames from there)
        autorun     - a flag for running clusterer automatically (on creation of object) 
    Note: The static method parse() can be used on any Robetta clusterer outfile, and returns a RobettaCluster
    object containing all cluster centers (RobettaClusterCenter objs) and a convergence values object (RobettaConvergence)
    """
    
    def __init__(self, command_file, decoy_file, sequence=None, outfile="data.cluster", \
                 log_file=None, executable="robetta_cluster_mod", path=None, autorun=True):
        self.command_file = command_file
        self.decoy_file   = decoy_file
        self.sequence     = sequence
        self.outfile      = outfile
        self.log_file     = log_file
        self.executable   = executable
        self.path         = path

        if self.path:
            os.chdir(self.path)
        if not self.sequence:
            self.sequence = self._get_sequence_from_file()
        if autorun:
            self.run()
    
    def run(self, ):
        # Create and populate robetta command file and run the clusterer
        self._create_command_file()
        self._execute()
        
        # Check for creation of clusterer output file
        if not os.path.isfile(self.outfile):
            raise Exception("Clusterer run failed, outfile '{0}' does not exist".format(self.outfile))

    @staticmethod
    def parse(result_file, centers=25):
    # Parse given Robetta clusterer result file. May return fewer than 'centers' (if file is short).
    # Returns a RobettaCluster object, containing a RosettaConvergence obj (convergence) and a list and dict of
    # RosettaClusterCenter objs (centers). See class RobettaCluster for details
        
        # NOTE: Cluster center (cc) pattern for HPF2 only
        header_pattern = r"THRESHOLD:\s+(?P<threshold>[0-9.]+)\s+TOP_CLUSTER_SIZE:\s+(?P<top_cluster_size>\d+)"   
        convergence_pattern = r"standard_thresholds:\s+size1=\s+(?P<size1>\d+)\s+threshold1=\s+(?P<t1>[0-9.]+)\s+size2=\s+(?P<s2>\d+)\s+threshold2=\s+(?P<t2>[0-9.]+)\s+total_decoys=\s+(?P<total_decoys>\d+)"
        cc_pattern = r"(?P<rank>\d+):\s+(?P<index>\d+),.+(?P<rosetta_id>[FS]_\d{4}_?\d{4})_?\s+(?P<size>\d+)\s+.*"
        stop_pattern = r"CLUSTER\sMEMBERS:"
        
        cluster_centers = list()

        with open(result_file) as handle:
            # Match first two lines against header and convergence values
            found = re.match(header_pattern, handle.next())
            if not found:
                raise Exception("Header not found in given Robetta clusterer outfile {0}".format(result_file))
            threshold, top_cluster_size = found.groups()

            found = re.match(convergence_pattern, handle.next())
            if not found:
                raise Exception("Convergence values not found in given Robetta clusterer outfile {0}".format(result_file))
            size1, radius1, size2, radius2, total_decoys = found.groups()
            convergence = RobettaConvergence(result_file, size1, radius1, total_decoys, size2=size2, radius2=radius2)

            # Match rest of file until 'CLUSTER MEMBERS:' with center pattern
            for line in handle:
                if re.match(stop_pattern, line):
                    break
                elif centers <= 0:
                    break
                found = re.match(cc_pattern, line)
                if found:
                    rank, index, rosetta_id, size = found.groups()
                    index = int(index)
                    cluster_centers.append( RobettaClusterCenter(rank, index, rosetta_id, size) )
                    centers -= 1
            cluster = RobettaCluster(result_file, threshold, top_cluster_size, convergence=convergence, centers=cluster_centers)
        return cluster
   

    def _get_sequence_from_file(self, ):
    # Parse sequence from given decoy file (in SILENT format). Sequence line should be first in file.
        handle = open(self.decoy_file)
        sequence = None
        for line in handle:
            found = re.match(r"SEQUENCE:\s*(?P<sequence>[A-Z]+)$", line)
            if found:
                sequence = found.group('sequence')
                break
        handle.close()
        if not sequence:
            raise Exception("Can not parse sequence from decoy file {0}".format(self.decoy_file))
        return sequence
    
    def _create_command_file(self, ):
    # Creates a robetta command file in CWD to be passed to robetta executable
    # Command file format:
    #OUTPUT_FILE <output>                - the file that clusterer outputs results to
    #TARGET <decoy_file> <sequence>      - the silent file containing decoys to cluster, and decoy sequence
        handle = open(self.command_file, 'w')
        handle.write("OUTPUT_FILE {0}\n".format(self.outfile))
        handle.write("TARGET {0} {1}\n".format(self.decoy_file, self.sequence))
        handle.close()

    def _execute(self, ):
        import subprocess
        cmd = [self.executable, self.command_file]
        if self.log_file:
            with open(self.log_file, 'w') as log_handle:
                ret = subprocess.check_call(cmd, stdout=log_handle, stderr=subprocess.STDOUT)
        else:
            ret = subprocess.check_call(cmd, stdout=None, stderr=None)
        return ret


class RobettaClustererParser(object):
    """ DEPRECATED """
    """ (for obvious reasons, I think) """

    HEADER = "header"
    RADII = "radii"
    RADIUS1 = "threshold1"
    RADIUS2 = "threshold2"
    SIZE1 = "size1"
    SIZE2 = "size2"
    SIZE_TOP = "size_top"
    TOTAL_DECOYS = "total_decoys"

    
    def parse(self,handle):
        from pyparsing import Group,Literal,lineEnd,empty
        from hpf.parsing import real,integer
        header = Group(Literal("THRESHOLD:")+
                       real.setResultsName(RobettaClustererParser.RADIUS1)+
                       Literal("TOP_CLUSTER_SIZE")+
                       real.setResultsName(RobettaClustererParser.SIZE_TOP)
                       ).setResultsName(RobettaClustererParser.HEADER)
        radii = Group(Literal("standard_thresholds:")+
                      Literal(RobettaClustererParser.SIZE1+"=")+
                      real.setResultsName(RobettaClustererParser.SIZE1)+
                      Literal(RobettaClustererParser.RADIUS1+"=")+
                      real.setResultsName(RobettaClustererParser.RADIUS1)+
                      
                      Literal(RobettaClustererParser.SIZE2+"=")+
                      real.setResultsName(RobettaClustererParser.SIZE2)+
                      Literal(RobettaClustererParser.RADIUS2+"=")+
                      real.setResultsName(RobettaClustererParser.RADIUS2)+
                      Literal(RobettaClustererParser.TOTAL_DECOYS+"=")+
                      integer.setResultsName(RobettaClustererParser.TOTAL_DECOYS)
                      ).setResultsName(RobettaClustererParser.RADII)
        
        expression = (header+lineEnd+
                      radii+lineEnd+
                      empty
                      )

        expression.parseFile(handle)
        assert False, "Save the data"


class RobettaClusterRun(Process):
    """ DEPRECATED """

    def __init__(self, control_file, log=None, bin="robetta_cluster_mod", **kwargs):
        super(RobettaClusterRun,self).__init__(**kwargs)
        self.control_file = control_file
        self.log = log
        self.bin = bin
        
    def _args(self):
    # Must pass shell=True to RosettaCluster (subprocess.Popen) to execute correctly
        cmd = "{0} {1}".format(self.bin, self.control_file)
        if self.log:
            cmd += " &> {0}".format(self.log)
        return cmd


class RobettaClusterCommandFile(ControlFile):
    """ DEPRECATED """
    
    def __init__(self, file, target, sequence, output="data.cluster", **kwargs):
        super(RobettaClusterCommandFile, self).__init__(file,**kwargs)
        self.target = os.path.abspath(target)
        self.output = os.path.abspath(output)
        self.sequence = sequence
    
    def _write(self, handle):
        print >>handle, "OUTPUT_FILE {0}".format(self.output)
        print >>handle, "TARGET {0} {1}".format(self.target, self.sequence)


def test():
    # Test and show usage of Robetta classes
    cmd_file = "na977.cluster.cmd"              # will create
    decoy_file  = 'na977.top100'                # must exist (silent file w/decoys)
    sequence    = "GVSCYGKPGGRLSGQLGESGEKGTNSSPGTRRFEQVTNSGRRGIFAKVEFVRVCVEYVYVVSYSLGFGVVLGVGLEYESDSILGSLINRGVCLL"
    cluster_out = "na977.top100.cluster.out"
    cluster_log = "na977.top100.log"
    
    print "Creating Clusterer obj and running clusterer"
    #RobettaClusterer (command_file, decoy_file, sequence=None, outfile='data.cluster', log_file=None, executable='robetta_cluster_mod', path=None, autorun=True)
    clusterer = RobettaClusterer(cmd_file, decoy_file, sequence=sequence, outfile=cluster_out, log_file=cluster_log)
    #clusterer = RobettaClusterer(cmd_file, decoy_file, path='/home/dpb3/mcm/newtest')

    print "Parsing clusterer output"
    cluster_results = RobettaClusterer.parse(clusterer.outfile)
    for cluster_center in cluster_results:
        print cluster_center

    #print "Extracting decoys as PDB files"
    #for cluster_center in cluster_results:
        # Run extractor (decoy_file, cluster_center.rosetta_id)

    print "Complete"

if __name__ == "__main__":
    test()



