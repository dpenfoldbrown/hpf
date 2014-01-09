
from hpf.hddb.db import url


from sqlalchemy.sql import func, and_
import sqlalchemy
import getopt
import sys
from numpy import average

import rpy2.robjects as robj

engine = sqlalchemy.create_engine(url+'pwinters')
metadata = sqlalchemy.MetaData(engine)
pdomainsccs = sqlalchemy.Table('p_domain_sccs',metadata, autoload=True)

hpf_engine = sqlalchemy.create_engine(url+'hpf')
hpf_metadata = sqlalchemy.MetaData(hpf_engine)
scop_des = sqlalchemy.Table('scop_des',hpf_metadata, autoload=True)

EUKARYOTIC_EXPERIMENTS = [804,886,827,826,825,888]
PROKARYOTIC_EXPERIMENTS = [889,840,809,814]
CUSTOM = [1176,]
NUM_SAMPLES = 10000
TINY_NUM = 0.0000000001


## dpb - I usually run (eg Berliner enrichments for RNA against whole-human):
# python sf_enrichment.py -x <hpf seqkey file> -s <exp>
# python sf_enrichment.py -x novel_hpf_seqids.txt -s 804
##

def usage():
        print "experiment=, sampling= (eukaryotic, prokaryotic)"

def main():
    print "welcome to main"
    experiment_keys = []
    sampling_experiments = "eukaryotic"
    filter_size = 0
    sequence_key_file = ""

    try:
        print "before getopt"
        opts, args = getopt.getopt(sys.argv[1:], "he:s:f:x:", ["help", "experiment=", "sampling=", "filter_size=", "sequence_key_file="])
        print "after getopt"
        print opts
        print args
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        print opt
        print arg
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-e","--experiment"):
            experiment_keys.append(arg)
        elif opt in ("-f","--filter_size"):
            filter_size = int(arg)
        elif opt in ("-x","--sequence_key_file"):
            sequence_key_file = arg
        elif opt in ("-s","--sampling"):
            if arg == "eukaryotic":
                sampling_experiments = EUKARYOTIC_EXPERIMENTS
            elif arg == "prokaryotic":
                sampling_experiments = PROKARYOTIC_EXPERIMENTS
            elif arg == "custom":
                sampling_experiments = CUSTOM
            else:
                sampling_experiments = map(int,arg.split(','))
        
        else:
            assert False, "unhandled option"


    #kdrew: query for superfamily prediction results
    sfp=[]
    if "" != sequence_key_file:
        seq_keys = []
        handle = open(sequence_key_file)
        for line in handle.readlines():
            seq_keys.append(line.rstrip())
        sfp.extend(get_superfamily_predictions_seqkey(seq_keys))
    sfp.extend(get_superfamily_predictions(experiment_keys))
    print "SFP enrichment set: {0}".format(len(sfp))

    #sfp_sampling = get_superfamily_predictions_sampling(EUKARYOTIC_EXPERIMENTS,NUM_SAMPLES,iterations=10)
    print "Sampling experiments: ", sampling_experiments
    sfp_sampling = [get_superfamily_predictions(sampling_experiments)]
    print "SFP sampling set: {0}".format(len(sfp_sampling))


    #kdrew: compile counts
    sfp_counts, sfp_psi_counts, sfp_fr_counts, sfp_denovo_counts = count_sf(sfp)
    sfp_counts = filterCountBySize(sfp_counts, filter_size)
    sfp_sampling_counts, sfp_sampling_psi_counts, sfp_sampling_fr_counts, sfp_sampling_denovo_counts = count_sf(sfp_sampling[0])
    sfp_sampling_counts = add_psuedo_count(sfp_sampling_counts, sfp_counts.keys())

    #kdrew: compile proportion
    sfp_proportions = calc_sf_proportion(sfp_counts)
    sfp_sampling_proportions = calc_sf_proportion(sfp_sampling_counts)

    #kdrew: compute enrichment
    sfp_enrichment = calc_enrichment(sfp_proportions, sfp_sampling_proportions)

    sfp_pvalue = calc_pvalue(sfp_counts, sfp_sampling_counts, bonferroni=False)
    sfp_pvalue_corrected = calc_pvalue(sfp_counts, sfp_sampling_counts, bonferroni=True)
    
    sorted_enrichment = sorted(sfp_pvalue.items(), key=lambda k: k[1])
    for i in sorted_enrichment:
        try:
            print i[0],";", sfp_enrichment[i[0]],";", sfp_counts[i[0]],";", sfp_sampling_counts[i[0]],";", scopDesc(i[0]),";", \
                    sfp_psi_counts[i[0]],";", sfp_fr_counts[i[0]],";", sfp_denovo_counts[i[0]],";", sfp_pvalue[i[0]],";", sfp_pvalue_corrected[i[0]]
        except KeyError:
            print i[0],";", sfp_enrichment[i[0]],";", sfp_counts[i[0]],";", 0,";", scopDesc(i[0]),";", \
                    sfp_psi_counts[i[0]],";", sfp_fr_counts[i[0]],";", sfp_denovo_counts[i[0]],";", sfp_pvalue[i[0]],";", sfp_pvalue_corrected[i[0]]
                
    
    #items = sfp_enrichment.items()
    #sorted_enrichment = sorted(items,key=lambda k: k[1])

    #j = len(sorted_enrichment)
    #for i in sorted_enrichment:
    #    try:
    #        sfp_sampling_counts[i[0]]
    #        print j,"    ", i[0],"    ", i[1],"    ", sfp_counts[i[0]],"    ", sfp_sampling_counts[i[0]],"    ", scopDesc(i[0]),"    ",sfp_psi_counts[i[0]],"    ",sfp_fr_counts[i[0]],"    ",sfp_denovo_counts[i[0]],"    ",sfp_pvalue[i[0]],"    ",sfp_pvalue_corrected[i[0]]
    #    except KeyError:
    #        print j,"    ", i[0],"    ",i[1],"    ", sfp_counts[i[0]],"    ", 0,"    ", scopDesc(i[0]),"    ",sfp_psi_counts[i[0]],"    ",sfp_fr_counts[i[0]],"    ",sfp_denovo_counts[i[0]],"    ",sfp_pvalue[i[0]],"    ",sfp_pvalue_corrected[i[0]]
    #
    #    j-=1

def add_psuedo_count(count_dict, keys, c=1):
    for key in keys:
        try:
            count_dict[key]+=c
        except KeyError:
            count_dict[key] = c

    return count_dict
        

#kdrew: only show superfamilies greater than or equal to 'size'
def filterCountBySize(unfiltered_dict, size):
    print "filtering by size %s" % size
    new_dict = dict()
    for key in unfiltered_dict.keys():
        if unfiltered_dict[key] >= size:
            new_dict[key] = unfiltered_dict[key]    

    return new_dict

def scopDesc(sccs):
    squery = sqlalchemy.select([scop_des.c.eng_desc], scop_des.c.sccs == sccs) 
    desc = squery.execute()
    return desc.fetchone()

#kdrew: convert a scop family sccs into its superfamily sccs
def family2superfamily(sccs):
    return '.'.join(sccs.strip().split('.')[:-1])

def family2fold(sccs):
    return '.'.join(sccs.strip().split('.')[:-2])

def get_superfamily_predictions_seqkey(sequence_keys,prob_min=0.8):
    pred_list = []
    squery = sqlalchemy.select([pdomainsccs.c.sccs,pdomainsccs.c.domain_type], sqlalchemy.and_(pdomainsccs.c.sccs != "NULL", pdomainsccs.c.parent_sequence_key.in_(sequence_keys), pdomainsccs.c.confidence >= prob_min) )
    preds = squery.execute()
    for pred in preds:
        print "sequence_key prediction %s" % pred
        for sccs in pred[0].split(','):
            pred_list.append((family2superfamily(sccs),pred[1]))

    return pred_list

def get_superfamily_predictions(experiments,prob_min=0.8):
    pred_list = []
    squery = sqlalchemy.select([pdomainsccs.c.sccs,pdomainsccs.c.domain_type], sqlalchemy.and_(pdomainsccs.c.sccs != "NULL", pdomainsccs.c.id.in_(experiments), pdomainsccs.c.confidence >= prob_min) )
    preds = squery.execute()
    for pred in preds:
        print "experiment prediction %s " % pred
        for sccs in pred[0].split(','):
            pred_list.append((family2superfamily(sccs),pred[1]))

    return pred_list
    
def get_superfamily_predictions_sampling(experiments,num_samples,iterations=1000):
    sample_list = []
    preds = get_superfamily_predictions(experiments)
    
    for i in xrange(1,iterations):
        print i
        sample_list.append(list(sample(preds, num_samples)))

    return sample_list

def calc_enrichment(proportion_dict1, proportion_dict2, log=True):
    enrichment_dict = dict()
    for key in proportion_dict1:
        try:
            enrichment_dict[key] = proportion_dict1[key] / proportion_dict2[key]
        except KeyError:
            #kdrew: this just means the sampling dict does not have the superfamily and therefore it is high enrichment
            enrichment_dict[key] = proportion_dict1[key] / TINY_NUM

        if log:
            import math
            enrichment_dict[key] = math.log(enrichment_dict[key])

    return enrichment_dict

def calc_pvalue(counts_dict1, counts_dict2, bonferroni=False):
    pvalue_dict = dict()
    #kdrew: find the total number of counts in the set (dict1) and total set+nonset (dict2)
    set_sum = 0 
    total_sum = 0
    bonferroni_count = 0
    for i in counts_dict1:
        set_sum += counts_dict1[i]
        bonferroni_count +=1
    for i in counts_dict2:
        total_sum += counts_dict2[i]

    print "set_sum: %s" % set_sum
    print "total_sum: %s" % total_sum
    print "bonferroni_count: %s" % bonferroni_count
    for key in counts_dict1:
        #kdrew: set up r matrix, assumes set (dict1) is a subset of total (dict2)
        sf_set = counts_dict1[key]
        sf_nonset = counts_dict2[key] - counts_dict1[key]
        #kdrew: error check, should never be but came across d.330.1 having this condition
        if sf_nonset < 0:
            continue
        nonsf_set = set_sum - counts_dict1[key]
        nonsf_nonset = total_sum - sf_nonset
        print "key: %s sf_set: %s sf_nonset: %s nonsf_set: %s nonsf_nonset: %s" % (key, sf_set, sf_nonset, nonsf_set, nonsf_nonset)
        v = robj.IntVector([sf_set, nonsf_set, sf_nonset, nonsf_nonset])
        m = robj.r['matrix'](v,nrow=2)
        ftest = robj.r['fisher.test']
        ft = ftest(m, alternative="greater")
        if bonferroni:
            pvalue_dict[key] = ft[0][0]*bonferroni_count
        else:
            pvalue_dict[key] = ft[0][0]

            

    return pvalue_dict

def calc_sf_proportion(count_dict):
    #kdrew: sum all domains in count_dict and normalize each entry
    total = 0
    proportion_dict = dict()
    for sf in count_dict.keys():
        total += count_dict[sf]

    print "total %s" % total

    for sf in count_dict.keys():
        proportion_dict[sf] = 1.0*count_dict[sf]/total

    return proportion_dict

def count_sf(pred_list):
    count_dict = dict()
    fr_dict=dict()
    psi_dict=dict()
    denovo_dict=dict()
    for pred in pred_list:
        try:
            count_dict[pred[0]] += 1
        except KeyError:
            count_dict[pred[0]] = 1
            fr_dict[pred[0]] = 0
            psi_dict[pred[0]] = 0
            denovo_dict[pred[0]] = 0

        if "fold_recognition" == pred[1]:
            fr_dict[pred[0]] += 1
        elif "psiblast" == pred[1]:
            psi_dict[pred[0]] += 1
        #else "unassigned" == pred[1] or "msa" == pred[1] or "pfam" == pred[1]:
        else: 
            denovo_dict[pred[0]] += 1

    return count_dict, psi_dict,fr_dict,denovo_dict


#kdrew: l is a list and n is number of samples
def sample(l,n):
    from numpy import array
    from numpy.random import shuffle

    x = array(l)
    shuffle(x)
    return x[:n]



if __name__ == "__main__":
    main()


