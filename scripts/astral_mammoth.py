#!/usr/bin/env python
# Mammoths all vs all astral domains.

import sys, os
import getopt
from hpf.pdb.mammoth import Mammoth,MammothCL,MammothAlignment
from hpf.utilities import consume
from hpf.processing import processor
from hpf.processing import MultiProcessor
from hpf.pdb.astral import AstralList
from script import paths
import MySQLdb

DEBUG = "debug"
SILENT = "silent"
UPLOAD = "noupload"
HELP = "help"
ASTRAL_DIR = "astral_dir"
opts = {DEBUG: False, SILENT: False, UPLOAD: True, ASTRAL_DIR:"/tmp/astral/"}
__hpf = None

def mammoth_set(id1):
	#kdrew: given id1 (the task) find all models that still need to be compared
	with _hpf(db="pdb") as cursor:
		query = 'select astral_id from astral95_1_75 as a where astral_id not in (select prediction from mammoth where experiment = "'+id1+'") and astral_id not in (select experiment from mammoth where prediction = "'+id1+'")'
		cursor.execute(query)
		compare2models = [t[0] for t in cursor.fetchall()]
	__hpf = None

	mm_list = []
	for id2 in compare2models:
		task = (id1,id2)
		#upload(mammoth(task))
		#mammoth(task)
		mm_list.append(mammoth(task))
		if len(mm_list) > 100:
			upload(mm_list)
			mm_list = []
	if len(mm_list) > 0:
		upload(mm_list)

def mammoth(task):
	id1,id2 = task
	print id1, id2
	al = AstralList()
	astral1 = al.retrieve_astral_file(id1)
	astral2 = al.retrieve_astral_file(id2)
	print astral1,astral2
	cl = MammothCL(astral1, astral2)
	m = Mammoth(cl)
	malign = m.run()
	malign.prediction = os.path.basename(astral1).split(".")[0]
	malign.experiment = os.path.basename(astral2).split(".")[0]
	return malign

def upload(mm_list):
	"""
	Uploads a mammoth alignment object
	@type mm: MammothAlignment
	"""
	print "uploading"
	db = _hpf(db="pdb")
	cursor = db.cursor()
	query = "insert into mammoth (prediction, experiment, ini_psi, ini_rms, end_psi, end_rms, zscore, evalue) values (%s,%s,%s,%s,%s,%s,%s,%s) on duplicate key update ini_psi=values(ini_psi),ini_rms=values(ini_rms),end_psi=values(end_psi),end_rms=values(end_rms),timestamp=NOW(),zscore=values(zscore),evalue=values(evalue)"

	upload_list = []
	for mm in mm_list:        
		upload_list.append((mm.prediction,mm.experiment,mm.ini_psi,mm.ini_rms,mm.end_psi,mm.end_rms,mm.zscore,mm.evalue))

	cursor.executemany(query,upload_list)
	cursor.close()
    
def retrieve(astral_id):
    al = AstralList()
    al.retrieve_astral_file(astral_id)

def tasks():
	global __hpf
	#tasks = []

	with _hpf(db="pdb") as cursor:
		query = "select astral_id from astral95_1_75"
		cursor.execute(query)
		models = [t[0] for t in cursor.fetchall()]
	__hpf = None
    
	print len(models)," models found"
	return models
	#models = list(paths.find(".*\.ent",dir=opts[ASTRAL_DIR]))

	#kdrew: tasks_list is a list of lists, each list goes to a job
	#tasks_list = []
	#for i,astral1 in enumerate(models):
	#	#kdrewa: task_set is a list of model pairs
	#	task_set = []
	#	for astral2 in models[i+1:]:
	#		task_set.append((astral1,astral2))
	#		print i, astral1, astral2
	#	tasks_list.append(task_set)
	#return tasks_list
    
def main():

    #pool = MultiProcessor(8,raise_errors=True,modulus=10)
    #for r in pool.run(retrieve,models):
    #    pass
    
    #models = models[:10]
    #for i,astral1 in enumerate(models):
    #    for astral2 in models[i+1:]:
    #        tasks.append((astral1,astral2))
    
    #tasks = tasks[1879500:]
    #tasks = tasks[200001:2000010]
    #pool = MultiProcessor(8, raise_errors=True)
    pool = processor(synchronous=False)
    #runtime().debug("Using processor",pool)
    pool.make_tasks(tasks)
    consume(pool.run(mammoth_set))
    #consume(pool.run(mammoth, tasks, upload))
    #for r in pool.run(mammoth, tasks, upload):
    #    pass
    print "FINISHED"

def _hpf(db=None):
    global __hpf
    if not __hpf:
        __hpf = MySQLdb.connect(db=db,read_default_file="~/.my.cnf")
    return __hpf

def _do(argv):
    try:
        opts, args = _options(argv)
    except:
        _usage()
        raise
    main(*args)

def _options(argv):
    _opts, args = getopt.getopt(argv, "?dxq", [DEBUG,SILENT,UPLOAD])
    global opts
    for o,a in _opts:
        if o in ('-?',HELP):
            _usage()
            sys.exit()
        if o in ('-d',DEBUG):
            opts[DEBUG] = True
        if o in ('-x', UPLOAD):
            opts[UPLOAD] = False
        if o in ('-q',SILENT):
            opts[SILENT] = True
        if o in ('',ASTRAL_DIR):
            opts[ASTRAL_DIR] = a
            assert os.path.exists(a)
    return opts, args

def _debug(*args):
    if opts[DEBUG]:
        _print(*args)

def _print(*args):
    if not opts[SILENT]:
        print " ".join([str(a) for a in args])

def _usage():
    print "$ template.py [-options] args"
    print "Template for quick scripts."
    print ""
    print "Options"


if __name__=="__main__":
	print "main"
	_do(sys.argv[1:])
