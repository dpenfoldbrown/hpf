#!/usr/bin/python

## make ALIGN file for folding homologues
## read argv[1] = psiblast file
## argv[2] = optional coords file (or '-')
## argv[3...] = list of .out files

import string
from os import popen,system,chdir,remove,getcwd
import sys




psiblast_file = sys.argv[1]
coords_file = sys.argv[2]
silent_files = sys.argv[3:]

#sys.stderr.write(string.join(silent_files)+'\n')




data = open(psiblast_file,'r')

query_sequence = {}

line = data.readline()
psiblast_output = 1

while line and line[:6] != 'Result':line = data.readline()
if not line:
    ## we're not in a psi-blast file
    psiblast_output = 0
    data.close()
    data = open(psiblast_file,'r')
    line = data.readline()

while line:
    if psiblast_output:
        assert line[:6] == 'Result'
        round_count = int(string.split(line)[3])

    align2query = {}
    full_seq = {}
    e_value = {}
    line = data.readline()
    while line and line[0] != '>':line = data.readline()
    while line and line[0] == '>':
        id = string.split(line)[0][1:]
        while line and line[:3] != ' Sc':line = data.readline()
        #repeat = 0 
        while line and line[:3] == ' Sc':
            #if repeat:print 'repeat',id
            #repeat = repeat +1
            e = string.split(line)[-1]
            if e[0] == 'e':e = '1'+e
            e = float(e)
            line = data.readline()
            qstart = 1000000
            hstart = 1000000
            qseq = ''
            hseq = ''
            while line and line[:5] != 'Query':line = data.readline()
            while line and line[:5] == 'Query':
                qstart = min(qstart,int(string.split(line)[1])-1)
                qstop = int(string.split(line)[3])-1
                qseq = qseq + string.split(line)[2]
                line = data.readline()
                line = data.readline()
                hstart = min(hstart,int(string.split(line)[1])-1)
                hstop = int(string.split(line)[3])-1
                hseq = hseq + string.split(line)[2]
                while line and line[:5] != 'Query':
                    if line[:3] == ' Sc' or \
                       line[0] == '>' or \
                       line[:6] == 'Result':break
                    line = data.readline()
            qs = string.join(string.split(qseq,'-'),'')
            ## testing query start
            #print qstop,qstart,len(qs),qs
            assert qstop-qstart+1==len(qs)
            assert len(hseq) == len(qseq)
            
            for i in range(len(qs)):
                pos = i + qstart 
                if pos not in query_sequence.keys():
                    query_sequence[pos] = qs[i]
                else:
                    assert query_sequence[pos] == qs[i]

            align = {}
            for i in range(len(hseq)):
                if hseq[i] == '-':continue
                pos = i - string.count(hseq[:i],'-')
                if qseq[i] == '-':
                    align[pos] = -1
                else:
                    qpos = i - string.count(qseq[:i],'-') + qstart
                    assert query_sequence[qpos] == qseq[i]
                    align[pos] = qpos
            
            hs = string.join(string.split(hseq,'-'),'')
            assert len(hs) == len(align.keys())
            
            unique_id = id+'_'+str(hstart)+'-'+str(hstop)
            align2query[unique_id] = align
            full_seq[unique_id] = hs
            e_value[unique_id] = e
    ks = query_sequence.keys()
    ks.sort()
    full_seq['QUERY'] = string.join(map(lambda x:query_sequence[x],ks),'')
    align2query['QUERY'] = {}
    for i in range(len(ks)):
        align2query['QUERY'][i] = ks[i]
        

    #sys.stderr.write(sys.argv[1]+' Round: '+`round_count`+\
    #                 ' Nseq: '+`len(e_value.keys())`+'\n')
    while line and line[:6] != 'Result':line = data.readline()
data.close()


file2id = {}
file_sequence = {}

for file in silent_files:
    silent_sequence = string.split(open(file,'r').readline())[1]
    file_sequence[file] = silent_sequence

    longest = 0
    for id in full_seq.keys():
        if string.count(silent_sequence,full_seq[id]) and len(full_seq[id])>longest:
            longest = len(full_seq[id])
            file2id [file] = id


    if not file2id.has_key( file ):
        sys.stderr.write('couldnt find sequence match for file '+file+\
                         '\nskipping!!!!!\n')

if coords_file != '-':
    line = open(coords_file,'r').readline()
    assert line[0] == '#'
    sequence = line[1:-1]
    file_sequence['NATIVE'] = sequence
    for id in full_seq.keys():
        if full_seq[id] == sequence:
            file2id ['NATIVE'] = id
            break

    if not file2id.has_key('NATIVE'):
        sys.stderr.write('couldnt find sequence match for coords file '+\
                         '\nskipping!!!!!\n')

    
    
used = []
for id in file2id.values():
    if id not in used:used.append(id)

## make a full alignment

gap_length = {} ## positions in blast query followed by a gap
ks = query_sequence.keys()
ks.sort()

gappy = {}
alignment = {}
for id in used:
    align = align2query[id]
    al = {} #alignment from query back
    for i in ks:al[i] = -1
    gp = {}
    prev = -1
    gap = []
    
    for pos in range(len(full_seq[id])):
        if align[pos] == -1:
            gap.append(pos)
        else:
            if gap:
                if gap_length.has_key(prev):
                    gap_length[prev] = max(len(gap),gap_length[prev])
                else:
                    gap_length[prev] = len(gap)
            gp[prev] = gap
            gap =[]
            al[ align[pos] ] = pos

            prev = align[pos]


    gappy[id] = gp
    alignment[id] = al

full_alignment = {}
for id in used:
    fa = ''
    al = alignment[id]
    gp = gappy[id]
    seq = full_seq[id]
    for pos in ks:
        if al[pos] == -1:
            fa = fa +'-'
        else:
            fa = fa + seq [ al[pos] ]
        if pos in gap_length.keys():
            gap = ''
            if pos in gp.keys():
                gap = string.join (map(lambda x: seq[x],gp[pos]),'')
            gap = gap + '-'*(gap_length[pos] - len(gap))
            fa = fa + gap
    full_alignment[id] = fa

## add terminal regions
L = len(full_alignment[file2id.values()[0]])

        
for file in file2id.keys():
    full_alignment [file] = full_alignment[file2id[file]]
    assert len(full_alignment[file]) == L
    
for file in file2id.keys():
    s1 = file_sequence[file]
    s2 = full_seq [ file2id [file] ]
    assert string.count(s1,s2)
    pos = string.find(s1,s2)
    L = L+pos
    ## add stuff at the N-terminus
    full_alignment[file] = file_sequence[file][:pos] + full_alignment[file]
    for file2 in file2id.keys():
        if file2 == file: continue 
        full_alignment[file2] = '-'*pos + full_alignment[file2]
        
    term = len(s1) - len(s2) - pos
    if term:
        L = L+term
        full_alignment[file] = full_alignment[file] + file_sequence[file][-1*term:]
        for file2 in file2id.keys():
            if file2 == file: continue 
            full_alignment[file2] =   full_alignment[file2] + '-'*term
        


for file in file2id.keys():
    assert len( full_alignment[file] ) == L
    print 'ALIGN',full_alignment [file],file    
    

   
