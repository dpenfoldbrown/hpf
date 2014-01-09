
import random
import cPickle as pickle
from Bio.PDB import *
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

# Debug flag. True prints, False does not.
DEBUG = False

#kdrew: this was taken from the methods section of Ora Schueler-Furman and David Baker (2003)
def inverse_distance_mean(dist_list):
    inverse_dist_list = map(lambda y: 1.0/y, dist_list)
    sum_idl = sum(inverse_dist_list)
    Ms = sum_idl/len(dist_list)

    if DEBUG: print "Inverse distance mean: {0}".format(Ms)
    if DEBUG: print "Length of given list: {0}".format(len(dist_list))
    return Ms

def get_average(list):
    sum = 0
    for el in list:
        sum += el
    return float(sum)/float(len(list))

def get_idm_percentile(target_idm, sample_set):
    greater = 0
    for sample_idm in sample_set:
        if sample_idm <= target_idm:
            greater += 1
    return float(greater) / float(len(sample_set)) * 100


class BioPDBStruct():
# Class to wrap opening a PDB file (by pdb id, eg 1a23) as a Bio.PDB object, and do structural
# operations on the PDB object (distance btw. residues, etc). Same functionality as below in
# the PymolPDBStruct class, but uses the cleaner/faster Bio.PDB (avoids opening Pymol)

    def __init__(self, pdb_id, chain, debug=False):
    # self.structure    - Bio.PDB pdb structure object
    # self.residues     - dict of Bio.PDB residue objects (triples of (heter_flag, sequence_number, insertion_code)),
    #                     of the form (res_seq_# => res_object)
        self.pdb_id = pdb_id
        self.chain  = chain
        self.DEBUG  = debug

        # Fetch and parse PDB file. Store in self.structure
        pdbl = PDBList()
        try:
            pdb_file = pdbl.retrieve_pdb_file(pdb_id)
        except:
            print "Retrieving PDB file {0} failed".format(pdb_id)
            raise
        parser = PDBParser()
        self.structure = parser.get_structure('local_struct', pdb_file)
        self.residues = self._get_clean_residues(self.structure, chain)

    def cluster_analysis(self, sites, sample_size=100, store_file=None, tag="", report=False):
    # sites         - a list of pdb atom numbers to consider in regards to 3d clustering
    # sample_size   - the number of random sample sets to run
    # store_file    - a pickle filename in which to store results. Will be appended to
    # tag           - a string containing any identifying information, used to store in in storefile
    #                 (if given) and print in report (if reporting)
    # report        - a boolean flag for printing statistical clustering analysis
        
        if sites == None:
            print "No sites provided for cluster analysis. Choosing sites at random"
            sites = self._get_random_sites()
        else:
            sites = self._cull_sites(sites)
        if len(sites) == 1:
            raise Exception("Only 1 site in sites list. Can not compute distances between 1 site")
        try:
            sites_distances = self.structure_distance_calc(sites)
        except Exception as e:
            print "Can not compute distance list for given sites"
            raise
        sites_idm = inverse_distance_mean(sites_distances)
        rand_idms = self.random_structure_distance_calc(len(sites), sample_size=sample_size)

        if store_file != None:
            self._store_cluster(tag, sites_idm, rand_idms, store_file)
        if report:
            self._report_cluster_stats(tag, sites_idm, rand_idms, sample_size)

    def structure_distance_calc(self, residues):
    # Returns a list of pairwise distances (in Angstroms) btw the given sites (ids)
        distance_list = []
        for i, residue1 in enumerate(residues):
            for residue2 in residues[(i+1):]:
                res1_atom = self._get_measurement_atom(residue1)
                res2_atom = self._get_measurement_atom(residue2)
                distance = res1_atom - res2_atom
                if distance == 0.0:
                    raise Exception("Distance between residues {0} and {1} returns {2}: Check pdb file".format(residue1, residue2, distance))
                distance_list.append(distance)
        return distance_list
    
    def random_structure_distance_calc(self, num_of_test_residues, sample_size=10):
    # Returns a list of IDMs of randomly selected sample sets of residues from the protein
        #if self.DEBUG: print "Residues from Bio.PDB structure: ", self.residues
        
        random_IDM_list = []
        for sample in range(0,sample_size):
            sample_ids = random.sample(self.residues.keys(), num_of_test_residues)
            if self.DEBUG: print "Randomly sampled residues: ", sample_ids
            try:
                sample_distance_list = self.structure_distance_calc(sample_ids)
            except Exception as e:
                print e
                print "Dropping this set from analysis. Moving to next sample set.."
                continue
                
            sample_IDM = inverse_distance_mean(sample_distance_list)
            random_IDM_list.append(sample_IDM)
        return random_IDM_list

    
    def _get_clean_residues(self, structure, chain):
    # Takes a Bio.PDB structure object.
    # Returns a dict of residue objects for a given chain WITHOUT water residues
    # Assumes structure model 0
        if not structure[0].has_id(chain):
            raise("Struct {0} does not have a chain '{1}'".format(structure.id, chain))
        res_dict = {}
        for res in structure[0][chain]:
            # If res is not a water, append to res_list (res.id format: triple, (hetero flag, seq #, ins code))
            if res.id[0] == 'W':
                continue
            res_dict[res.id[1]] = res
        return res_dict

    def _get_measurement_atom(self, residue):
    # Takes a Bio.PDB residue object (triple)
    # Returns a Bio.PDB atom object representing an atom in the residue to provide a 3d-space measurement center.
        if self.residues[residue].get_resname() == 'GLY':
            atom = 'CA'
        else:
            atom = 'CB'
        try:
            return self.residues[residue][atom]
        except KeyError as e:
            print "Residue {0} does not contain requested atom {1}".format(residue, atom)
            raise

    def _get_random_sites(self, ):
    # Return a set of sites from the structure randomly.
        num_residues = len(self.residues)
        if num_residues < 100:
            num_sites = 5
        elif num_residues < 200:
            num_sites = 10
        else:
            num_sites = 15
        return random.sample(self.residues.keys(), num_sites)
    
    def _cull_sites(self, sites):
    # Remove all sites from sites list not in the object's structure (self.residues)
    # return a new list of culled sites
        culled_sites = []
        for residue in sites:
            if residue not in self.residues:
                continue
            culled_sites.append(residue)
        return culled_sites

    def _store_cluster(self, tag, idm, rand_idms, store_file='cluster_records.default.pkl'):
    # Check and open (w append) given store file, and write cluster info to it
    # Written format is a quintuple, of (given tag info, pdb id, pdb chain, idm of sites, randomly sampled idms)
        pkl_handle = open(store_file, 'ab')
        pickle.dump((tag, self.pdb_id, self.chain, idm, rand_idms), pkl_handle)
        pkl_handle.close()

    def _report_cluster_stats(self, tag, idm, rand_idms, sample_size):
    # Print cluster statistical information to screen
        print "Structure {0}{1}, {2}".format(self.pdb_id, self.chain, tag)
        print "   IDM for given sites: {0}".format(idm)
        print "   Average randomly sampled IDM: {0}, from {1} samples".format(get_average(rand_idms), sample_size)
        print "   Sites' IDM greater than {0} percent of randomly sampled sets' IDMs".format(get_idm_percentile(idm, rand_idms))




class PymolPDBStruct():
# Class to wrap opening a PDB-formatted structure file as a Pymol object, and to do
# structural operations on the Pymol object.

    def __init__(self, structure, debug=False):
    # structure - a structure filesting in PDB file format. (Filestring - contents of file as a string).
        import pymol as pm
        self.pm = pm
        self.pm.finish_launching()
        self.pm.cmd.read_pdbstr(structure, "struct")
        
        # Use pymol to get all residue ids (kdrew). 
        self.pm.stored.resids = []
        self.pm.cmd.iterate("struct and name ca", "stored.resids.append(resi)")
        self.DEBUG = debug

        # Do some magic (convert to set => removes duplicates) to remove duplicates from list
        # (Because some PDB files return duplicates, for fun).
        # Convert resid list to ints, then cast as a set (to remove duplicates), sort (for debugging only), and assign.
        self.pm.stored.resids = sorted(list(set(map(int, self.pm.stored.resids))))
        if self.DEBUG: print "Residues taken from structure via pymol: ", self.pm.stored.resids
    
    def clear(self,):
    # Wipes everything from the pymol object (which, for some reason, persists in memory like crazy).
        self.pm.cmd.delete('all')

    #kdrew: returns list of Ms for randomly selected sets of residues
    def random_structure_distance_calc(self, num_of_test_residues, sample_size=10):
        
        if self.DEBUG: print "Resids from pymol: ", self.pm.stored.resids
        random_Ms_list = []
        for sample in range(0,sample_size):
            sample_ids = random.sample(self.pm.stored.resids, num_of_test_residues)
            if self.DEBUG: print "Randomly sampled residues: ", sample_ids
            try:
                rand_d_list = self.structure_distance_calc(sample_ids)
            except Exception as e:
                print e
                print "Dropping this set from analysis. Moving to next sample set.."
                continue
                
            rand_Ms = inverse_distance_mean(rand_d_list)
            random_Ms_list.append(rand_Ms)
    
        return random_Ms_list

    #kdrew: not sure of the residue indexing openbabel uses, i.e. starts at 0 through num_of_residues?, does not seem to follow pdb indexing
    #kdrew: maybe pymol was the better choice for this indexing reason
    def structure_distance_calc(self, ids):
        
        distance_list = []
        #kdrew: iterate through all pairs of residue ids
        for i, resid1 in enumerate(ids):
            for resid2 in ids[(i+1):]:
                self.pm.cmd.select("res1_c", "struct and (name cb or (name ca and resn GLY)) and resi %s" % resid1) 
                self.pm.cmd.select("res2_c", "struct and (name cb or (name ca and resn GLY)) and resi %s" % resid2)
                
                #DEBUG
                #self.pm.cmd.iterate("res1_c", "print chain,resi,name")
                #self.pm.cmd.iterate("res2_c", "print chain,resi,name")
    
                distance = self.pm.cmd.distance("dist_res1_res2", "res1_c", "res2_c")
                if distance == 0.0:
                    raise Exception("Distance between residues {0} and {1} returns {2}: Check pdb file".format(resid1, resid2, distance))
                #if self.DEBUG: print "Residues: %s , %s Distance: %f" % (resid1, resid2, distance)
                distance_list.append(distance)
    
        #if self.DEBUG: print "List of distances: ", distance_list
        return distance_list

