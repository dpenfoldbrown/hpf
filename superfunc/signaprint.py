#!/usr/bin/env python

# A class and parsing file for creating signaprints (signatures ofr mouse sequences
# based on MGI, GO terms, IPRs, Pfams, and Gene names)

import argparse
import cPickle as pickle
from hpf.interpro.interpro import InterproResultFile, InterproResult

class Signaprint():
    
    def __init__(self, sequence_id, sequence_length=None, mgi=None, gene_sym=None, iprs=None, pfams=None, go_mfs=None, go_bps=None, go_ccs=None):
    # iprs, pfams, and all go variables are lists
        self.sequence_id = sequence_id
        self.sequence_length = sequence_length
        self.mgi = mgi
        self.gene_syms = gene_sym
        self.iprs = iprs
        self.pfams = pfams
        self.go_mfs = go_mfs
        self.go_bps = go_bps
        self.go_ccs = go_ccs
        self.go_leaves = None

    def __repr__(self, ):
        rep = "<Signaprint: {0}".format(self.sequence_id)
        if self.mgi:
            rep += ", {0}".format(self.mgi)
        if self.iprs:
            rep += " IPRs: {0}".format(len(self.iprs))
        if self.pfams:
            rep += " Pfams: {0}".format(len(self.pfams))
        if self.go_mfs and self.go_bps and self.go_ccs:
            rep += " GO terms: {0}".format(len(self.go_mfs) + len(self.go_bps) + len(self.go_ccs))
        return rep

    def get_all_GO(self, ):
    # Returns a list of ALL unique go terms associated with this print
    # Returns empty list if no go terms present
        all_go = list()
        if self.go_mfs != None:
            all_go += self.go_mfs
        if self.go_bps != None:
            all_go += self.go_bps
        if self.go_ccs != None:
            all_go += self.go_ccs
        all_go = list(set(all_go))
        return all_go


def create_signaprints():
    # Call all file parse methods to populate Signaprint object dict (seq_key => SignaPrint)
    signaprints = init_dict()

    populate_mgi(signaprints)
    populate_seqlength(signaprints)
    populate_iprs(signaprints, '/home/dpb/superfunc/data/interpro/mouse1171_interpro_all.out')
    populate_pfams(signaprints, '/home/dpb/superfunc/data/pfam/mouse1171_pfam_all.out')
    
    populate_GO(signaprints, '/home/dpb/superfunc/data/GO/mouseGoMF.txt', "MF")
    populate_GO(signaprints, '/home/dpb/superfunc/data/GO/mouseGoBP.txt', "BP")
    populate_GO(signaprints, '/home/dpb/superfunc/data/GO/mouseGoCC.txt', "CC")
    populate_GO_leaves(signaprints)

    print "Creating signaprint dict complete. Entries: {0}".format(len(signaprints))
    
    outfile = 'mouse_signaprints.pkl'
    handle = open(outfile, 'wb')
    pickle.dump(signaprints, handle)
    handle.close()


def populate_GO_leaves(sdict):
# Takes list of all GOs (all types) and runs the leaf-finder, setting Signaprint.go_lowest to the list of localleaf GO terms
    from hpf.superfunc.go import local_leaf_finder
    have_goleaves = 0
    for seq in sdict.keys():
        all_goterms = sdict[seq].get_all_GO()
        if not all_goterms: continue
        sdict[seq].go_leaves = local_leaf_finder(all_goterms)
        have_goleaves += 1
    print "Parsed local GO leaves for {0} seqs".format(have_goleaves)

def populate_GO(sdict, filename, type):
    handle = open(filename)
    # If first line is key, drop it
    handle.next()
    for line in handle:
        fields = line.rstrip().split('\t')
        seq_id = int(fields[0])
        gene_sym = fields[1]
        mgi_id = fields[4]
        is_not = int(fields[7])
        go_id = fields[9]
        
        if is_not == 1:
            go_id = "NOT_" + go_id
        # Add GO symbol to correct list
        if type == "MF":
            if not sdict[seq_id].go_mfs:
                sdict[seq_id].go_mfs = [go_id]
            elif go_id not in sdict[seq_id].go_mfs:
                sdict[seq_id].go_mfs.append(go_id)
        elif type == "BP":
            if not sdict[seq_id].go_bps:
                sdict[seq_id].go_bps = [go_id]
            elif go_id not in sdict[seq_id].go_bps:
                sdict[seq_id].go_bps.append(go_id)
        elif type == "CC":
            if not sdict[seq_id].go_ccs:
                sdict[seq_id].go_ccs = [go_id]
            elif go_id not in sdict[seq_id].go_ccs:
                sdict[seq_id].go_ccs.append(go_id)
        # Add gene_name
        if not sdict[seq_id].gene_syms:
            sdict[seq_id].gene_syms = [gene_sym]
        elif gene_sym not in sdict[seq_id].gene_syms:
            sdict[seq_id].gene_syms.append(gene_sym)

def init_dict():
# Creates a returns a dict of the form sequence id => Signaprint obj. Initially, signaprint object is empty
    from hpf.hddb.db import Session, Protein
    from sqlalchemy import distinct
    session = Session()
    mouse_seqs = session.query(distinct(Protein.sequence_key)).filter_by(experiment_key=1171)
    signaprint_dict = dict()
    for (seq,) in mouse_seqs:
        signaprint_dict[seq] = Signaprint(sequence_id=seq)
    print "Dict of {0} sequence => Signaprints created".format(len(signaprint_dict))
    return signaprint_dict

def populate_mgi(sdict):
    from hpf.hddb.db import Session, SequenceAc
    session = Session()
    mgi_count = 0
    for seq in sdict.keys():
        acs = session.query(SequenceAc).join(SequenceAc.protein).filter_by(sequence_key=seq, experiment_key=1171).all()
        mgi_ids = set()
        for ac in acs:
            if ac.ac2 == "None":
                continue
            mgi_ids.add(ac.ac2)
        if len(mgi_ids) == 0:
            sdict[seq].mgi = "None"
        elif len(mgi_ids) == 1:
            sdict[seq].mgi = list(mgi_ids)[0]
            mgi_count += 1
        else:
            sdict[seq].mgi = list(mgi_ids)
            mgi_count += 1
    print "{0} sequences assigned an MGI ID".format(mgi_count)

def populate_seqlength(sdict):
    from hpf.hddb.db import Session, Sequence
    session = Session()
    for seq_id in sdict.keys():
        seq_obj = session.query(Sequence).get(seq_id)
        sdict[seq_id].sequence_length = len(seq_obj.sequence)
    print "Getting sequence length for all sequences complete"

def populate_iprs(sdict, ipr_resultfile):
    interpro_results = InterproResultFile(filename=ipr_resultfile)
    have_iprs = 0
    for seq in interpro_results.results_byquery:
        iprs = set()
        for result in interpro_results[seq]:
            if result.ipr != "NULL":
                iprs.add(result.ipr)
        sdict[int(seq)].iprs = list(iprs)
        if len(iprs) > 0: have_iprs += 1
    del(interpro_results)
    print "{0} sequences assigned list of IPRs".format(have_iprs)

def populate_pfams(sdict, pfam_resultfile):
    pfam_results = InterproResultFile(filename=pfam_resultfile)
    have_pfams = 0
    for seq in pfam_results.results_byquery:
        pfams = set()
        for result in pfam_results[seq]:
            pfams.add(result.result)
        sdict[int(seq)].pfams = list(pfams)
        if len(pfams) > 0: have_pfams += 1
    del(pfam_results)
    print "{0} sequences assigned list of Pfams".format(have_pfams)

def main():
    if args.create:
        create_signaprints()
        print "Creating signaprint dict complete"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Signaprint creation and handling')
    parser.add_argument('-c', action="store_true", default=False, dest='create', help='Create a Signaprint dict by calling create() function')
    args = parser.parse_args()
    main()

