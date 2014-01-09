'''
Created on Feb 24, 2010

@author: patrick
'''
from hpf.runtime import runtime
from hpf.align import Alignment, TemporaryAlignmentFile, AlignmentFactory
from hpf.seq import TemporaryRecordFile
from itertools import imap
from tempfile import NamedTemporaryFile
import subprocess

from hpf.mapper import MapperInterface

class CulledColumnMapper(MapperInterface):
    """
    Yields a (culled_column, original_column) mapper.  Useful for mapping CodeML
    results (run on a culled alignment) to the original family.
    """
    
    def __init__(self, alignment, culled_columns):
        """
        @type alignment: Bio.Align.Generic.Alignment, hpf.hddb.db.Alignment
        @type alignment: iterable of hpf.hddb.db.CulledColumn or int
        """
        super(CulledColumnMapper,self).__init__()
        from Bio.Align.Generic import Alignment
        self.alignment = alignment if isinstance(alignment,Alignment) else alignment.alignment
        self.culled_columns = [c.column if hasattr(c,"column") else c for c in culled_columns]
        
    def mapping(self):
        culled = set([c for c in self.culled_columns])
        i = 0
        for column in xrange(self.alignment.get_alignment_length()):
            if column in culled:
                continue
            yield (i,column)
            i=i+1
            

class SeedAlignment(Alignment):
    """
    Represents an alignment where the final record is considered to be the
    target added to an unchanged seed alignment.
    """
    
    def __init__(self, *args, **kwargs):
        Alignment.__init__(self,*args, **kwargs)
        self._mapping = {}

    def __len__(self):
        return len(self._records)

    def _copy(self, alignment):
        Alignment._copy(self, alignment)
        self._targets = [record for record in self._records if not record.id.startswith("_seed_")]
        self._rename()

    def _rename(self):
        for record in self._records:
            record.id = record.id.replace("_seed_","")

    def map(self, target_num=0):
        """
        Map residue numbers in the target record to columns from the original
        seed alignment. Caches mapping for future calls.
        ie.
        >seed1
        A-BCD
        >seed2
        A-BCD
        >target
        -ZBCD
        returns [None,1,2,3]
        
        @return: list, The original alignment column that the target mapped to.
            len(return)==len(target.seq), beginning with column 0.
        """
        if self._mapping.has_key(target_num):
            return self._mapping[target_num]
        _mapping = []
        seed_col = 0
        is_gapped = lambda x: x=='-'
        for col in xrange(self.get_alignment_length()):
            aa = self.get_column(col)
            #print "TARGETNUM",target_num,"TARGETS",self._targets
            assert len(self._targets)>target_num
            target_is_gapped = is_gapped(self._targets[target_num].seq._data[col])
            # If the seed alignment is completely gapped at this column,
            # then this represents an insertion by the target... meaning the
            # target doesn't align on this column.
            seed_is_gapped = all(imap(is_gapped, aa[0:len(self)-len(self._targets)]))
            #print map(is_gapped, aa[0:len(self)-len(self._targets)]),len(self)-len(self._targets)
            # if the target record is not gapped, determine what seed column
            # it maps to.  If the seed alignment is gapped here, then it
            # doesn't align and returns None
            if not target_is_gapped:
                #print seed_is_gapped,seed_col,len(_mapping)
                _mapping.append(None if seed_is_gapped else seed_col)
            # If the seed alignment is not gapped, we increment the column
            # tracker
            if not seed_is_gapped:
                seed_col+=1
        # The mapping should only be as long as the target record
        assert all([len(str(target.seq).replace("-","")) == len(_mapping) for target in self._targets])
        self._mapping[target_num] = _mapping
        return _mapping

    def seed_mapping(self):
        is_gapped = lambda x: x=='-'
        original_column = 0
        for col in xrange(self.get_alignment_length()):
            aa = self.get_column(col)
            # If the seed alignment is completely gapped at this column,
            # then this represents an insertion by the target... meaning the
            # target doesn't align on this column.
            seed_is_gapped = all(imap(is_gapped, aa[0:len(self)-len(self._targets)]))
            if not seed_is_gapped:
                yield original_column
                original_column+=1
            else:
                yield None
    
    def seed(self):
        return self._records[0:len(self)-len(self._targets)]
    
    def targets(self):
        return self._targets

class SeedAlignmentFactory(object):
    
    def __init__(self, _class=SeedAlignment):
        self._class = _class

    def create(self, seed, targets, format="fasta", **kwargs):
        """
        Performs a default Mafft seed alignment using the seed and target
        @return: Seed alignment object
        """
        output = NamedTemporaryFile(**kwargs)
        with TemporaryAlignmentFile([seed], format=format, **kwargs) as seed_file:
            with TemporaryRecordFile(targets, format=format, **kwargs) as target_file:
                cmd = "mafft-linsi --quiet --seed %s %s > %s" % (seed_file,target_file, output.name)#"test.txt"
                runtime().debug(cmd)
                subprocess.check_call(cmd,shell=True)
        #print open(output.name).read()
        with open(output.name) as handle:
            alignment = AlignmentFactory(self._class).read(handle,format)
        return alignment
    
