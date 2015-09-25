#!/usr/bin/python

import string
from glob import glob
from os.path import exists
from math import floor,log10,log,exp
from operator import add
from os import chdir,system,popen
import sys
from whrandom import random
from math import sqrt

#print '\n\nTHIS SCRIPT IS UNDER CONSTRUCTION --- SORRY'

MAX_RMS = 15

if len(sys.argv)<2:
    print '\n'+'-'*75
    print 'usage:',sys.argv[0],' {-gs} {-e} {-ss <ss-prefix>} {-f <frag-file>} <contact-file1> {<contact-file2> ... }\n\n\n'
    print '-gs gives grey-scale output (default is rainbow colors)'
    print '-e specifies .eps format (better for combining multiple plots)'
    print '-ss allows you to pass a prefix (ss-prefix) to which the program'
    print '  will add .rdb, .phd, .psipred_ss2, and .jufo_ss in looking for'
    print '  secondary structure information'
    print '-f will read the ss-info from the designated frag-file\n'+'-'*75+'\n\n'
    assert 0==1

######### parse args
arguments = sys.argv[1:]

if '-gs' in arguments:
    pos = arguments.index('-gs')
    del arguments[pos]
    GREY_SCALE = 1
else:
    GREY_SCALE = 0


if '-e' in arguments:
    pos = arguments.index('-e')
    del arguments[pos]
    format = 'eps'
else:
    format = 'ps'

if '-ss' in arguments:
    pos = arguments.index('-ss')
    ss_base = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]
else:
    ss_base = ''

if '-f' in arguments:
    pos = arguments.index('-f')
    frag_file = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]
else:
    frag_file = ''
    
    
data_files = arguments

SKIP_NEARBY = 5

def Fig(xpos,ypos,marked_color,color,width,marked,out,height=0.0):
    if height==0.0:height=width
    ypos = 10000-ypos
    if marked:
        out.write(string.join(map(str,['2 2 0 1',marked_color,color,
                                       '50 0 20 0.000 0 0 -1 0 0 5\n',
                                       '\t',xpos,ypos,
                                       width+xpos,ypos,
                                       width+xpos,height+ypos,
                                       xpos,height+ypos,
                                       xpos,ypos]))+'\n')
    else:
        out.write(string.join(map(str,['2 2 0 0',color,color,
                                       '50 0 20 0.000 0 0 -1 0 0 5\n',
                                       '\t',xpos,ypos,
                                       width+xpos,ypos,
                                       width+xpos,height+ypos,
                                       xpos,height+ypos,
                                       xpos,ypos]))+'\n')
def Line(xpos,ypos,width,out):
    ypos = 10000-ypos
    out.write(string.join(map(str,['2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n',
                                   xpos,ypos,xpos+width,ypos]))+'\n')
    return

def RMS100(rmsd,length):
    return rmsd / (1 + 0.5 * (log ( length/100.0))) 

def Text(xpos,ypos,label,out):
    ypos = 10000-ypos
    out.write(string.join(map(str,['4 0 0 50 0 0 12 0.0000 4 195 135',
                                   xpos,ypos,label+'\\001']))+'\n')


def Color(fraction,LOG):
    if LOG:
        EXPONENT = log10(2)
    else:
        EXPONENT = 1.0
        
    if fraction < 0.00001:
        score = fraction
    else:
        score = exp( EXPONENT * log(fraction))

    color = 32 + int(floor(2*16*16*(score/1.00001)))
    return color


def Color_RMSD(rmsd):
    if rmsd > MAX_RMS:
        fraction = 0.0
    else:
        fraction = (float(MAX_RMS) - rmsd) / MAX_RMS
    return Color(fraction, 0) ## linear interpolation between 0 and 15

def Read_phd(phd_file):
    lines = map(string.split,open(phd_file,'r').readlines())
    e = []
    h = []
    l = []
    seq= '' 
    for line in lines:
        if len(line) == 2 and line[0] == 'AA' and line[1][0] == '|':
            seq = seq + line[1][1:-1]
        elif len(line) == 1 and line[0][:5] == 'prH-|':
            h = h + map(lambda x:float(x)/10,list(line[0][5:]))
        elif len(line) == 1 and line[0][:5] == 'prE-|':
            e = e + map(lambda x:float(x)/10,list(line[0][5:]))
        elif len(line) == 1 and line[0][:5] == 'prL-|':
            l = l + map(lambda x:float(x)/10,list(line[0][5:]))
    L = len(seq)
    if len(h) != L or len(e) != L or len(l) != L:
        sys.stderr.write('WARNING: error reading phd file: %s\n'%phd_file)
        return [[],[],[],'']
    return [e,h,l,seq]

def Read_fragments(fragment_file):
    sys.stderr.write('Reading fragment file: %s\n'%fragment_file)
    data = open(fragment_file,'r')
    line = data.readline()
    prev = (-1,-1)
    ss_count = {}
    while line:
        if len(line)>20 and line[-5] == 'F' and line[-10] == 'P':
            position = int(line[-9:-6])
            fragment = int(line[-4:-1])
            if fragment <= 25:
                if (position,fragment) != prev: ## first in fragment
                    prev = (position,fragment)
                    count = 0
                    if fragment == 1: ## initialize counts
                        if not position%10:
                            sys.stderr.write('read frag file: position %d\n'%position)
                        
                        for i in range(9):
                            fpos = position+i-1 ## numbering starts at 0
                            if fpos not in ss_count.keys():
                                ss_count[fpos] = {'H':0,'E':0,'L':0}
                    
                pos = position + count - 1## numbering starts at 0
                count = count+1

                ss = line[16]
                ss_count[pos][ss] = ss_count[pos][ss] + 1
        line = data.readline()
    data.close()
    
    l = []
    e = []
    h = []
    for i in range(len(ss_count.keys())):
        total = reduce(add,ss_count[i].values())
        if total>0:
            e.append( float(ss_count[i]['E'])/total)
            h.append( float(ss_count[i]['H'])/total)
            l.append( float(ss_count[i]['L'])/total)

    return e,h,l

def Read_jufo(jufo_file):
    lines = map(string.split,open(jufo_file,'r').readlines())
    l = []
    e = []
    h = []
    seq = ''
    for line in lines:
        if len(line)== 6:
            seq = seq + line[1]
            l.append(float(line[3]))
            h.append(float(line[4]))
            e.append(float(line[5]))
    return [e,h,l,seq]
    
def Read_jones(jones_file):
    lines = map(string.split,open(jones_file,'r').readlines())
    l = []
    e = []
    h = []
    seq = ''
    for line in lines:
        if len(line)== 6:
            seq = seq + line[1]
            l.append(float(line[3]))
            h.append(float(line[4]))
            e.append(float(line[5]))
    return [e,h,l,seq]
    
def Read_rdb(rdb_file):
    lines = open(rdb_file,'r').readlines()
    while lines and lines[0][0] == '#':del lines[0]
    del lines[:2]
    lines = map(string.split,lines)
    seq = ''
    e = []
    h = []
    l = []
    for line in lines:
        if len(line) != 5:continue
        seq = seq + line[1]
        e.append(float(line[2]))
        h.append(float(line[3]))
        l.append(float(line[4]))
    return [e,h,l,seq]

def Read_ss_info(ss_base):
    ss_info = {}

    rdb_file = ss_base+'.rdb'
    if exists(rdb_file):
        print 'reading rdb file:',rdb_file
        ss_info['rdb'] = Read_rdb(rdb_file)
    else:
        print 'couldnt find rdb file:',rdb_file
        
    jones_file = ss_base+'.psipred_ss2'
    if exists(jones_file):
        print 'reading jones file:',jones_file
        ss_info['jon'] = Read_jones(jones_file)
    else:
        new_jones_file = string.join(string.split(jones_file,'/')[:-1],'/')+'/psipred_ss2'
        if exists(new_jones_file):
            print 'WARNING: couldnt find',jones_file
            print 'WARNING: Using',new_jones_file,'instead'
            ss_info['jon'] = Read_jones(new_jones_file)
        else:
            print 'couldnt find jones file:',jones_file

    jufo_file = ss_base+'.jufo_ss'
    if exists(jufo_file):
        print 'reading jufo file:',jufo_file
        ss_info['juf'] = Read_jufo(jufo_file)
    else:
        print 'couldnt find jufo file:',jufo_file

    phd_file = ss_base+'.phd'
    if exists(phd_file):
        print 'reading phd file:',phd_file
        ss_info['phd'] = Read_phd(phd_file)
    else:
        print 'couldnt find phd file:',phd_file

    return ss_info

def Color_line(values,pos1,pos2,width,height,out):
    for i in range(len(values)):
        color = Color(values[i],0) ## not logarithmic
        if 1 or values[i]>0:
            Fig(pos1+width*i, pos2, color,color,width,0,out,height)
    return

        
def Contact_plot(ps_file,contact_file,info,format,LOG=1):
    ss_color = {'H':5,'E':6,'L':7}
    MARKED_COLOR = 0


    if ss_base:
        ss_info = Read_ss_info(ss_base)
    else:
        ss_info = {}
    
    if frag_file:
        frag_ss = Read_fragments(frag_file)
    else:
        frag_ss = {}

    lines = map(string.split,open(contact_file,'r').readlines())

    native_contacts = []
    native_ss = ''
    contact_fraction = {}
    maxsub_fraction = {}
    decoy_ss = [ [], [], [] ] ############## e, h, l
    ssl = ['E','H','L']
    subcluster_rmsd = {}
    subcluster_threshold = {}
    DISPLAY_SUBCLUSTERING = 0
    for line in lines:
        if line[0] == 'NS':
            native_ss = line[1]
        elif line[0] == 'NC':
            native_contacts.append( (int(line[1]),int(line[2])) )
        elif line[0] == 'DS' or line[0] == 'SS':
            pos = int(line[1])
            assert pos == len(decoy_ss[0])
            for i in range(3):
                ss = line[2*i+2][0]
                decoy_ss[ ssl.index(ss) ].append(float(line[2*i+3]))
        elif line[0] == 'DC' or line[0] == 'CC':
            assert 0<= float(line[3]) <= 1
            contact_fraction[(int(line[1]),int(line[2]))] = float(line[3])
        elif line[0] == 'MS':
            assert 0<= float(line[3]) <= 1
            maxsub_fraction[(int(line[1]),int(line[2]))] = float(line[3])
        elif line[0] == 'SR':
            subcluster_rmsd [(int(line[1]),int(line[2]))] = float(line[3])
            DISPLAY_SUBCLUSTERING = 1
        elif line[0] == 'ST':
            subcluster_threshold [(int(line[1]),int(line[2]))] = float(line[3])
            DISPLAY_SUBCLUSTERING = 1

    L = len(decoy_ss[0])

    ## add the 0.0 data:
    st = {}
    sr = {}
    cf = {}
    ms = {}
    marked = {}
    for i in range(L):
        for j in range(L):
            cf[(i,j)] = 0.0
            ms[(i,j)] = 0.0
            marked[(i,j)] = 0
            sr[(i,j)] = -1
            st[(i,j)] = -1
    
    for pos in contact_fraction.keys():
        assert pos[0]<L and pos[1]<L
        cf[pos] = contact_fraction[pos]
        
    for pos in maxsub_fraction.keys():
        assert pos[0]<L and pos[1]<L
        ms[pos] = maxsub_fraction[pos]

    for pos in subcluster_rmsd.keys():
        sr[pos] = subcluster_rmsd[pos]

    for pos in subcluster_threshold.keys():
        st[pos] = subcluster_threshold[pos]

    contact_fraction = cf
    maxsub_fraction = ms
    subcluster_rmsd = sr
    subcluster_threshold = st

    for pos in native_contacts:marked[pos] = 1
    
    if native_ss:
        assert len(native_ss) == L
    else:
        native_ss = 'L'*L
        

    fig_file = '/tmp/phil_junk'+str(random())+'.fig'

    sys.stderr.write('Making %s\n'%ps_file)
    out = open(fig_file,'w')
##     out.write(string.join(open('/users/pbradley/header.fig','r').readlines(),''))
    out.write('#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nSingle\n-2\n1200 2\n')

    ###########################  define new colors
    l = map(str,range(10))+['a','b','c','d','e','f']



    counter = 32
    if GREY_SCALE:
        for i in range(256):
            red = i/2
            green = i/2
            blue = i/2
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1

        for i in range(256):
            red = 128+i/2
            green = 128+i/2
            blue = 128+i/2
            
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1
    else:
        for i in range(256):
            if LOG:
                green = min(200,2*i)
                blue = min(255,510-2*i)
                red = 0
            else:
                green = 255- i/2
                red = 255-i/2
                blue = 255
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1

        for i in range(256):
            if LOG:
                red = min(255,2*i)
                green = min(200,510-2*i)
                blue = 0
            else:
                green = 128-i/2
                red = 128-i/2
                blue = 255
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1

    
    width = (12500/(2*L + 5))

    offset = width*(L+5)
    
    to_fig = []
    for i in range(L):
        for j in range(i+SKIP_NEARBY,L):
            k = (i,j)
            
            color = Color(contact_fraction[k],LOG)
            
            pos1 = i*width
            pos2 = j*width


            ## contacts below the diagonal
            if marked[k]:
                to_fig.append([pos2,pos1,MARKED_COLOR,color,width,1,out])
            else:
                Fig(pos2,pos1,color,color,width,0,out)


            ## maxsub above the diagonal
            color = Color(maxsub_fraction[k],LOG)
            
            Fig(pos1,pos2,color,color,width,0,out)

            if (k[1] == L-1 and not k[0]%10):
                Fig(pos1,pos2+width,0,0,width,0,out)
                if not k[0]%50:
                    Fig(pos1,pos2+2*width,0,0,width,0,out)

            ## subcluster_rmsd above diagonal, offset
            if subcluster_rmsd [k] != -1:

                color = Color_RMSD ( subcluster_rmsd[k])
                Fig (pos1+offset,pos2,color,color,width,0,out)

            if subcluster_threshold [k] != -1:
                color = Color_RMSD (subcluster_threshold[k])
                Fig (pos2+offset,pos1,color,color,width,0,out)




    for tf in to_fig:
        Fig(tf[0],tf[1],tf[2],tf[3],tf[4],tf[5],tf[6])


    text_list = [[0,(L+2)*width,info]]

    if 1:
        ## make secondary structure plot ####################33

        #calculate deviations of decoys from predictions
        dev = [ [0.5]*L, [0.5]*L ] # 0=E,1=H
        ks =ss_info.keys()
        for type in ks:
            if (len(ss_info[type][0]) !=L):
                sys.stderr.write('WARNING: prediction/sequence mismatch: %s\n'%type);
                del ss_info[type]
        if ss_info:
            for pos in range(L):
                for ss in range(2):
                    pred = 0.0
                    for type in ss_info.keys():
                        pred = pred + ss_info[type][ss][pos]
                    pred = pred/len(ss_info.keys())
                    dev[ss] [pos] = max(0.0,min(1.0, 0.5 + (decoy_ss[ss][pos]-pred)))

        ## smooth dev
        for i in range(2):
            for pos in range(L):
                if pos==0:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos+1])/2
                elif pos==L-1:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos-1])/2
                else:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos-1]+dev[i][pos+1])/3
        
        assert len(decoy_ss[0]) == len(native_ss)
        native = [[],[],[]]
        for i in range(L):
            if native_ss[i]=='E':
                native[0].append(1.0)
                native[1].append(0.0)
                native[2].append(0.0)
            elif native_ss[i]=='H':
                native[0].append(0.0)
                native[1].append(1.0)
                native[2].append(0.0)
            elif native_ss[i]=='L':
                native[0].append(0.0)
                native[1].append(0.0)
                native[2].append(1.0)

        if frag_ss:
            if len(frag_ss[0]) != L:
                sys.stderr.write('WARNING: fragment file length mismatch\n')
                frag_ss = {}

        pos2 = L*width+600

        text_down = -125
        text_over = -500
        line_up = 3
        DW=200
        PW=150
        for ss in range(2): ## 0=strand, 1=helix
            ss_name = 'EH'[ss]
            ## decoys:
            text_list.append([text_over,pos2+text_down,'dec: '+ss_name])
            Color_line(decoy_ss[ss],0,pos2,width,DW,out)
            Line(0,pos2+line_up,width*L,out)
            pos2 = pos2+PW+10

            ## fragments
            if frag_ss:
                text_list.append([text_over,pos2+text_down,'frg: '+ss_name])
                Color_line(frag_ss[ss],0,pos2,width,PW,out)
                Line(0,pos2+line_up,width*L,out)
                pos2 = pos2+PW+10
                

            ## predictions
            for type in ss_info.keys():
                if len(ss_info[type][ss]) == L:
                    text_list.append([text_over,pos2+text_down,type+': '+ss_name])
                    Color_line(ss_info[type][ss],0,pos2,width,PW,out)
                    pos2 = pos2+PW+10
            
            ## native
            Line(0,pos2-PW-10+line_up,width*L,out)
            text_list.append([text_over,pos2+text_down,'nat: '+ss_name])
            Color_line(native[ss],0,pos2,width,PW,out)
            pos2 = pos2+DW+50
        Color_line(dev[0],0,pos2,width,DW,out)
        text_list.append([text_over,pos2+text_down,'dev: E'])
        pos2 = pos2+DW
        Color_line(dev[1],0,pos2,width,DW,out)
        text_list.append([text_over,pos2+text_down,'dev: H'])
        
        ## consensus decoy ss
        pred_ss = ''
        for pos in range(L):
            if decoy_ss[0][pos] > decoy_ss[1][pos] and decoy_ss[0][pos] > decoy_ss[2][pos]:
                pred_ss = pred_ss + 'E'
            elif decoy_ss[1][pos] > decoy_ss[2][pos]:
                pred_ss = pred_ss + 'H'
            else:
                pred_ss = pred_ss + 'L'

        assert len(pred_ss) == L

        for i in range(len(pred_ss)):
            ss = pred_ss[i]

            if DISPLAY_SUBCLUSTERING:
                Fig(i*width+offset,i*width,ss_color[ss],ss_color[ss],width,0,out)

            Fig(i*width,i*width,ss_color[ss],ss_color[ss],width,0,out)
            Fig(i*width,(L+1)*width,ss_color[ss],ss_color[ss],width,0,out)
            Fig(-1*width,i*width,ss_color[ss],ss_color[ss],width,0,out)
            Fig(i*width,-1*width,ss_color[native_ss[i]],ss_color[native_ss[i]],
                width,0,out)
            Fig((L+1)*width,i*width,ss_color[native_ss[i]],
                ss_color[native_ss[i]],width,0,out)
            

    ## make key #######################3
    show = [1000,500,100,50,10,5,1]
    labels = {}
    for s in show:
        cf = float(s)/1000
        labels[Color(cf,LOG)] = str(s)


    w2 = 5000/512
    for i in range(512):
        color = 32+i

        pos1 = -150
        pos2 = i*w2

        for j in range(10):
            Fig(pos1+j*w2,pos2,color,color,w2,0,out)

        if color in labels.keys():
            for j in range(10):
                Fig(pos1-1-j*w2,pos2,color,color,w2,0,out)

            text_list.append([pos1-500,pos2,labels[color]])

    ## make key for contacts and maxsub #######################3
    show = [1000,500,100,50,10,5,1]
    labels = {}
    for s in show:
        cf = float(s)/1000
        labels[Color(cf,LOG)] = str(s)


    w2 = 5000/512
    for i in range(512):
        color = 32+i

        pos1 = -150
        pos2 = i*w2

        for j in range(10):
            Fig(pos1+j*w2,pos2,color,color,w2,0,out)

        if color in labels.keys():
            for j in range(10):
                Fig(pos1-1-j*w2,pos2,color,color,w2,0,out)

            text_list.append([pos1-500,pos2,labels[color]])

    
    if DISPLAY_SUBCLUSTERING: ## make key for rmsd100 #######################3

        lengths = [25,50,100,150]
        rmsds = range(1,15)
        labels  = {}
        for r in rmsds:
            for l in lengths:
                if RMS100(r,l) >= MAX_RMS:
                    continue
                labels [ Color_RMSD (RMS100(r,l)) ] = [r,l] 

        for l in lengths:
            pos1 = offset + 10*width - 10*w2
            pos2 = L*width + 2*width

            up = lengths.index(l)
            text_list.append([pos1,pos2+10*w2+200*up,str(l)])



        w2 = 5000/512
        for i in range(512):
            color = 32+i

            pos1 = offset + 10*width + (511-i)*w2
            pos2 = L*width + 2*width

            for j in range(10):
                Fig(pos1,pos2-j*w2,color,color,w2,0,out)

            if color in labels.keys():
                for j in range(10):
                    Fig(pos1,pos2+1+j*w2,color,color,w2,0,out)

                up = lengths.index(labels[color][1])
                text_list.append([pos1,pos2+10*w2+200*up,str(labels[color][0])])


    ## write all the text
    for tl in text_list:
        Text(tl[0],tl[1],tl[2],out)

    out.close()


    system('fig2dev -L '+format+' -l 0 '+fig_file+' '+ps_file)
    system('rm '+fig_file)
    
    return






    


for file in data_files:
    Contact_plot(file+'.'+format,file,file,format,1)
