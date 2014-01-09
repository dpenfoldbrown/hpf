#!/usr/bin/env python
import MySQLdb

db = MySQLdb.connect(host="mysql", user="pwinters", passwd="bonneaulab", db="pwinters")
cursor = db.cursor()

def ex(stmt):
    print stmt
    cursor.execute(stmt)

def drop(table):
    try:
        ex("drop table %s" % table)
    except:
        pass

drop("domain_sccs")
ex("create table domain_sccs like hpf.domain_sccs")
ex("insert ignore into domain_sccs select * from hpf.domain_sccs")

drop("experiment")
ex("create table experiment like hpf.experiment")
ex("insert ignore into experiment select * from hpf.experiment")

drop("protein")
ex("create table protein like bddb.protein")
ex("insert ignore into protein select p.* from experiment e join bddb.protein p on e.id=p.experiment_key")

drop("domain")
ex("create table domain like ddbCommon.domain")
ex("insert ignore into domain select d.* from protein p join ddbCommon.domain d on p.sequence_key=d.parent_sequence_key")

drop("mcmData")
ex("create table mcmData like bddb.mcmData")
ex("insert ignore into mcmData select m.* from domain d join bddb.mcmData m on d.domain_sequence_key=m.sequence_key")

drop("filesystemOutfile")
ex("create table filesystemOutfile like bddb.filesystemOutfile")
ex("insert ignore into filesystemOutfile select f.* from domain d join bddb.filesystemOutfile f on d.outfile_key=f.id")
ex("insert ignore into filesystemOutfile select f.* from domain d join bddb.filesystemOutfile f on d.domain_sequence_key=f.parent_sequence_key")
ex("insert ignore into filesystemOutfile select f.* from domain d join bddb.filesystemOutfile f on d.domain_sequence_key=f.sequence_key")

drop("sequence")
ex("create table sequence like ddbCommon.sequence")
for table in ["ddbCommon.ddb_sequence", "ddbCommon.sequence"]:
    ex("insert ignore into sequence select s.id,s.sha1,s.sequence from protein p join %s s on p.sequence_key=s.id" %table)
    ex("insert ignore into sequence select s.id,s.sha1,s.sequence from domain d join %s s on d.domain_sequence_key=s.id"%table)
    ex("insert ignore into sequence select s.id,s.sha1,s.sequence from domain d join %s s on d.parent_sequence_key=s.id"%table)
    ex("insert ignore into sequence select s.id,s.sha1,s.sequence from filesystemOutfile f join %s s on f.sequence_key=s.id"%table)

drop("sequenceAc")
ex("create table sequenceAc like ddbCommon.sequenceAc")
for table in ["ddbCommon.ddb_sequenceAc","ddbCommon.sequenceAc"]:
    ex("insert ignore into sequenceAc select a.* from sequence s join %s a on s.id=a.sequence_key" % table)

drop("nrTaxName")
ex("create table nrTaxName like ddbCommon.nrTaxName")
ex("insert ignore into nrTaxName select n.* from sequenceAc a join ddbCommon.nrTaxName n on a.taxonomy_id=n.taxonomy_id")

drop("nrTaxNode")
ex("create table nrTaxNode like ddbCommon.nrTaxNode")
ex("insert ignore into nrTaxNode select * from ddbCommon.nrTaxNode")

drop("go")
ex("""
CREATE TABLE `go` (
  `sequence_key` int(10) unsigned NOT NULL DEFAULT '0',
  `acc` varchar(10) DEFAULT NULL,
  `code` varchar(8) NOT NULL,
  KEY `sequence_key` (`sequence_key`),
  KEY `acc` (`acc`) );""")
ex("""
insert ignore into go SELECT sm.id as sequence_key,t.acc,e.code FROM ddbCommon.mygo_association AS ass
  INNER JOIN ddbCommon.mygo_gene_product AS gp ON ass.gene_product_id = gp.id
  INNER JOIN ddbCommon.mygo_gene_product_seq AS gps ON gp.id = gps.gene_product_id
  INNER JOIN ddbCommon.sequenceMeta as sm ON gps.seq_id = sm.mygo
  INNER JOIN ddbCommon.mygo_term AS t ON term_id = t.id
  INNER JOIN ddbCommon.mygo_evidence AS e ON association_id = ass.id
  INNER JOIN protein as p ON p.sequence_key = sm.id;""")


import subprocess
print "mysqldump -h mysql -u pwinters -pbonneaulab pwinters | gzip > hpf.sql.gz"
tables = "domain_sccs experiment protein domain mcmData filesystemOutfile sequence sequenceAc nrTaxName nrTaxNode"
retcode = subprocess.call("mysqldump -h mysql -u pwinters -pbonneaulab pwinters %s | gzip > /work1/wintersp/hpf.sql.gz" % tables, shell=True)
if retcode != 0:
    raise Exception("Subprocess error",retcode)
