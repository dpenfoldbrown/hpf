
#kdrew: go_hpf_utils.R is for methods which query the mygo/hpf database

source("sql_setup.R")

#kdrew: method for retreiving acc, name, term_type ... from term table in mygo database
get_goinfo <- function(id)
{
	query = paste("SELECT * from term where id = ", id);
	as.list(mygo_qry(query));
}



#kdrew: TODO change functionFrom to 'go_from' 
getFunction <- function(seq = NULL, acc=FALSE, name=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA','NR'), functionFrom, create_table=FALSE, seq_table=NULL, allmfs=FALSE)
{
	if(allmfs)
	{
		return(getFunction_allmfs(functionFrom=functionFrom, create_table=create_table));
	}
	else
	{
		return(getFunction_generic(seq=seq, acc=acc, name=name, code=code, functionFrom=functionFrom, create_table=create_table, seq_table=seq_table));
	}
}

#drop table if exists mfs; create temporary table mfs (select DISTINCT t.id as mf_id  from seq2go as s2g, sequence_class as sc, mygo.gene_product_seq as gps, mygo.term as t  where s2g.gene_product_id = gps.gene_product_id and gps.seq_id = sc.seq_id   and t.id = s2g.term1_id and t.term_type = 'molecular_function'   and s2g.code IN ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA')  and sc.seq_id = 40761);

#kdrew: this creates a table which lists ids for all molecular functions in the GO database (not just the ones annotated to the testSet)
getFunction_allmfs <- function(functionFrom="mygo", create_table=TRUE)
{
	query = "";
	
	if(create_table)
	{
		query_drop = paste("drop table if exists mfs");
		print(query_drop);
		hpf_qry(query_drop);
		
		query = paste("create temporary table mfs (");
	}
	
	query = paste(query, "select DISTINCT t.id as mf_id ");
	query = paste(query, " from ",functionFrom,".term as t ", sep="");
	query = paste(query, " where t.term_type = 'molecular_function'  )"); 
	
	print(query_drop);
	print(query);
	
	hpf_qry(query_drop);
	return(hpf_qry(query));
}

#kdrew: generic function for querying go annotations for a specific benchmark
getFunction_generic<- function(seq = NULL, acc=FALSE, name=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'), functionFrom, create_table = FALSE, seq_table=NULL)
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = "";
	
	if(create_table)
	{
		query_drop = paste("drop table if exists mfs");
		print(query_drop);
		hpf_qry(query_drop);
		
		query = paste("create temporary table mfs (");
	}
	
	query = paste(query, " select DISTINCT d.domain_sequence_key as seq_id, gp.term1_id as mf_id, t.acc as mf_acc ");
	if(name)
	{
		query = paste(query, ", t2.name");
	}
	query = paste(query, " from hpf.domain as d, ",functionFrom," as f, ", sep="");
	query = paste(query, " mygo.term as t, mygo.term as t2, mygo.graph_path as gp ");
	query = paste(query, " where d.parent_sequence_key = f.sequence_key ");
	query = paste(query, " and t.acc = f.acc and gp.term2_id = t.id and gp.term1_id = t2.id");
	query = paste(query, " and f.code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and t.term_type = 'molecular_function' and t.is_obsolete = 0");
	if(!is.null(seq))
	{
		query = paste(query, " and d.domain_sequence_key = ", seq );
	}
	if(!is.null(seq_table))
	{
		query = paste(query, " and d.domain_sequence_key in (select seq_id from ",seq_table, ")" );
	}
	
	if(create_table)
	{
		query = paste(query, " )");
	}
	
	print(query);
	funcs = hpf_qry(query);
	
	if(create_table)
	{
		alter_query = paste("alter table mfs add index(seq_id)");
		hpf_qry(alter_query);
		alter_query = paste("alter table mfs add index(mf_id)");
		hpf_qry(alter_query);
	}
	return(funcs);
}

getProcess <- function(seq, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA','NR'), processFrom = "goSet", create_table=FALSE, negative=FALSE)
{
	if("gold_standard" == processFrom)
	{ 
		return(getProcessGoldStandard(seq=seq, acc=acc, code=code));
	}
	if("scopfold" == processFrom)
	{ 
		return(getProcessScopfold(seq=seq, acc=acc, code=code));
	}
	if("goSet" == processFrom)
	{
		return(getProcessGoSet(seq=seq, acc=acc, code=code));
	}
	else
	{
		return(getProcess_generic(seq=seq, acc=acc, code=code, processFrom=processFrom, create_table=create_table, negative=negative));
	}
}

#kdrew: doesn't return acc
getProcessScopfold <- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT gp.term1_id as bp_id ");
	query = paste(query, " from scopfold.go_terms as f, mygo.graph_path as gp ");
	query = paste(query, " where gp.term2_id = f.term_id ");
	query = paste(query, " and f.evidence_code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and f.term_type = 'biological_process' ");
	if(!is.null(seq))
	{
		query = paste(query, " and f.sequence_key = ", seq );
	}
	
	print(query);
	procs = hpf_qry(query);
	return(procs);
	
}
#kdrew: doesn't return acc
getProcessGoldStandard <- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT gp.term1_id as bp_id ");
	query = paste(query, " from gold_standard.domain as d, gold_standard.function as f, mygo.graph_path as gp ");
	query = paste(query, " where d.sequence_key = f.sequence_key and gp.term2_id = f.term_id ");
	query = paste(query, " and f.evidence_code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and f.term_type = 'biological_process' ");
	if(!is.null(seq))
	{
		query = paste(query, " and d.domain_sequence_key = ", seq );
	}
	
	print(query);
	procs = hpf_qry(query);
	return(procs);
	
}
#kdrew: this function returns (or creates a table) with a list of the biological process terms associated with the given sequence including parent nodes
#processFrom is where to get the process terms, create_table if set to TRUE will create  temporary sql table, 
#negative if TRUE will return all process terms not associated with the given sequence
getProcess_generic <- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'), processFrom, create_table=FALSE, negative=FALSE)
{
	query = "";
	
	if(create_table)
	{
		query_drop = paste("drop table if exists neg_procs");
		print(query_drop);
		hpf_qry(query_drop);
		query_drop = paste("drop table if exists procs");
		print(query_drop);
		hpf_qry(query_drop);
		
		query = paste("create temporary table procs (");
	}
	
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste(query, "select DISTINCT gp.term1_id as bp_id ");s
	query = paste(query, " from hpf.domain as d, ",processFrom," as f,", sep=""); 
	query = paste(query, " mygo.term as t, mygo.graph_path as gp ");
	query = paste(query, " where t.acc = f.acc and d.parent_sequence_key = f.sequence_key and gp.term2_id = t.id");
	query = paste(query, " and f.code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and t.term_type = 'biological_process' and t.is_obsolete=0 ");
	if(!is.null(seq))
	{
		query = paste(query, " and d.domain_sequence_key = ", seq );
	}
	
	if(create_table)
	{
		query = paste(query, " )");
	}
	
	print(query);
	procs = hpf_qry(query);
	
	if(create_table)
	{
		#alter_query = paste("alter table procs add index(seq_id)");
		#hpf_qry(alter_query);
		alter_query = paste("alter table procs add index(bp_id)");
		hpf_qry(alter_query);
	}
	if(negative && create_table)
	{
		neg_query = paste("create table neg_procs"); 
		neg_query = paste(neg_query, "select id as bp_id from mygo.term where term_type = 'biological_process' and  is_obsolete = 0 and id not in (select bp_id from procs)");
		print(neg_query);
		neg_procs = hpf_qry(neg_query);
		
		alter_query = paste("alter table neg_procs add index(bp_id)");
		hpf_qry(alter_query);
		
		return(neg_procs);
	}
	
	return(procs);
	
}

getProcessGoSet <- function(seq, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT t.id as bp_id");
	if(TRUE == acc)
	{
		query = paste(query, ", t.acc as bp_acc ");
	}
	query = paste(query, " from seq2go as s2g, sequence_class as sc, mygo.gene_product_seq as gps, mygo.term as t");
	query = paste(query, " where s2g.gene_product_id = gps.gene_product_id and gps.seq_id = sc.seq_id ");
	query = paste(query, " and t.id = s2g.term1_id and t.term_type = 'biological_process' ");
	query = paste(query, " and s2g.code IN (", go_evidence, ")"); 
	query = paste(query, " and sc.seq_id = ", seq);
	procs = hpf_qry(query);
	return(procs);
}
getLocalization <- function(seq, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA','NR'), localizationFrom = "goSet", create_table=FALSE, negative=FALSE)
{
	if("gold_standard" == localizationFrom)
	{ 
		return(getLocalizationGoldStandard(seq=seq, acc=acc, code=code));
	}
	if("scopfold" == localizationFrom)
	{ 
		return(getLocalizationScopfold(seq=seq, acc=acc, code=code));
	}
	if("goSet" == localizationFrom)
	{
		return(getLocalizationGoSet(seq=seq, acc=acc, code=code));
	}
	else
	{
		return(getLocalization_generic(seq=seq, acc=acc, code=code, localizationFrom=localizationFrom, create_table=create_table, negative=negative));
	}
}
#kdrew: doesn't return acc
getLocalizationScopfold<- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT gp.term1_id as cc_id ");
	query = paste(query, " from scopfold.go_terms as f, mygo.graph_path as gp ");
	query = paste(query, " where gp.term2_id = f.term_id ");
	query = paste(query, " and f.evidence_code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and f.term_type = 'cellular_component' ");
	if(!is.null(seq))
	{
		query = paste(query, " and f.sequence_key = ", seq );
	}
	
	print(query);
	locs = hpf_qry(query);
	return(locs);
	
}
#kdrew: doesn't return acc
getLocalizationGoldStandard <- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT gp.term1_id as cc_id ");
	query = paste(query, " from gold_standard.domain as d, gold_standard.function as f, mygo.graph_path as gp ");
	query = paste(query, " where d.sequence_key = f.sequence_key and gp.term2_id = f.term_id ");
	query = paste(query, " and f.evidence_code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and f.term_type = 'cellular_component' ");
	if(!is.null(seq))
	{
		query = paste(query, " and d.domain_sequence_key = ", seq );
	}
	
	print(query);
	locs = hpf_qry(query);
	return(locs);
	
}

getLocalization_generic<- function(seq = NULL, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'), localizationFrom, create_table=FALSE, negative=FALSE)
{
	query = "";
	if(create_table)
	{
		query_drop = paste("drop table if exists neg_locs");
		print(query_drop);
		hpf_qry(query_drop);
		query_drop = paste("drop table if exists locs");
		print(query_drop);
		hpf_qry(query_drop);
		
		query = paste("create temporary table locs (");
	}
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste(query, "select DISTINCT gp.term1_id as cc_id ");
	query = paste(query, " from hpf.domain as d, ",localizationFrom," as f,",sep=""); 
	query = paste(query, " mygo.term as t, mygo.graph_path as gp ");
	query = paste(query, " where t.acc = f.acc and d.parent_sequence_key = f.sequence_key and gp.term2_id = t.id ");
	query = paste(query, " and f.code in (", go_evidence, ")");
	query = paste(query, " and gp.term1_id != 1 and t.term_type = 'cellular_component' and is_obsolete=0 ");
	if(!is.null(seq))
	{
		query = paste(query, " and d.domain_sequence_key = ", seq );
	}
	if(create_table)
	{
		query = paste(query, " )");
	}
	
	print(query);
	locs = hpf_qry(query);
	if(create_table)
	{
		#alter_query = paste("alter table locs add index(seq_id)");
		#hpf_qry(alter_query);
		alter_query = paste("alter table locs add index(cc_id)");
		hpf_qry(alter_query);
	}
	
	if(negative && create_table)
	{
		neg_query = paste("create table neg_locs"); 
		neg_query = paste(neg_query, "select id as cc_id from mygo.term where term_type = 'cellular_component' and  is_obsolete = 0 and id not in (select cc_id from locs)");
		print(neg_query);
		neg_locs = hpf_qry(neg_query);
		
		alter_query = paste("alter table neg_locs add index(cc_id)");
		hpf_qry(alter_query);
		
		return(neg_locs);
	}
	
	return(locs);
	
}
getLocalizationGoSet <- function(seq, acc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT t.id as cc_id ");
	if(TRUE == acc)
	{
		query = paste(query, ", t.acc as cc_acc ");
	}
	query = paste(query, " from seq2go as s2g, sequence_class as sc, mygo.gene_product_seq as gps, mygo.term as t");
	query = paste(query, " where s2g.gene_product_id = gps.gene_product_id and gps.seq_id = sc.seq_id ");
	query = paste(query, " and t.id = s2g.term1_id and t.term_type = 'cellular_component' ");
	query = paste(query, " and s2g.code IN (", go_evidence, ")"); 
	query = paste(query, " and sc.seq_id = ", seq);
	procs = hpf_qry(query);
	return(procs);
}

getStructure <- function(seq, sf=FALSE, structureFrom = "goSet", sf_data_table="mcm_data")
{
	if("gold_standard" == structureFrom)
	{ 
		return(getStructureGoldStandard(seq=seq, sf=sf));
	}
	if("scopfold" == structureFrom)
	{ 
		return(getStructureScopfold(seq=seq, sf=sf));
	}
	if("goSet" == structureFrom)
	{
		return(getStructureGoSet(seq=seq, sf=sf));
	}
	else
	{
		return(getStructure_generic(seq=seq, sf=sf, structureFrom=structureFrom, sf_data_table=sf_data_table));
	}
}

getStructureScopfold <- function(seq, sf=FALSE)
{
	query = paste("select sf_id from scopfold.mcm_data as mds where mds.sequence_key =", seq);
	print(query);
	sfs = hpf_qry(query);
	return(sfs);
}
getStructure_generic<- function(seq, sf=FALSE, structureFrom, sf_data_table="mcm_data")
{
	query = paste("select sf.id as sf_id from ",structureFrom,".",sf_data_table," as md, hpf.superfamily as sf ",sep="");
	query = paste(query, " where md.sccs = sf.superfamily and md.domain_sequence_key =", seq);
	print(query);
	sfs = hpf_qry(query);
	return(sfs);
}
getStructureGoldStandard <- function(seq, sf=FALSE)
{
	query = paste("select sf_id from gold_standard.mcm_data_scaled as mds where mds.domain_sequence_key =", seq);
	print(query);
	sfs = hpf_qry(query);
	return(sfs);
}

#kdrew: for a given sequence return its superfamily id
getStructureGoSet <- function(seq, sf=FALSE)
{
	query = paste("select gbt.sf_id as sf_id ");
	if(TRUE == sf)
	{
		query = paste(query, ", SUBSTRING_INDEX(gbt.sccs_id,'.',3) as superfamily ");
	}
	query = paste(query, "from go_blast_trim as gbt ");
	query = paste(query, "where gbt.seq_id = ", seq);
	sfs = hpf_qry(query);
	return(sfs);
}

getStructureProb <- function(seq, struct, structureProbFrom = "goSet")
{
	if("gold_standard" == structureProbFrom)
	{ 
		getStructureProbGoldStandard(seq=seq, acc=acc, code=code);
	}
	if("goSet" == structureProbFrom)
	{
		getStructureProbGoSet(seq=seq, acc=acc, code=code);
	}
}
getStructureProbGoSet <- function(seq, struct)
{
	query = paste("select ????");
	
}

getSeqsWithGO <- function(class=-1, mf=TRUE, bp=FALSE, cc=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA'), completed=TRUE)
{
	go_type = NULL;
	if(mf)
	{ 	
		go_type = c(go_type, "molecular_function")
	}
	if(bp)
	{ 	
		go_type = c(go_type, "biological_process")
	}
	if(cc)
	{ 	
		go_type = c(go_type, "cellular_component")
	}
	type_collapse = paste('\'',go_type,'\'', collapse=',', sep='');
	
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');
	
	query = paste("select DISTINCT sc.seq_id ");
	query = paste(query, " from seq2go as s2g, sequence_class as sc, mygo.gene_product_seq as gps, mygo.term as t");
	query = paste(query, " where s2g.gene_product_id = gps.gene_product_id and gps.seq_id = sc.seq_id ");
	query = paste(query, " and t.id = s2g.term1_id ");
	
	query = paste(query, " and t.name IN (", type_collapse, ")");
	query = paste(query, " and s2g.code IN (", go_evidence, ")"); 
	if(-1 != class)
	{
		query = paste(query, " and sc.class = ",class);
	}
	if(completed)
	{
		query = paste(query, " and sc.prob_computed = 1");
	}
	
	print(query);
	seqs = hpf_qry(query);
	return(seqs);
}

#kdrew: return seq ids in a given class, complement returns all seq ids not in given class
subsetSeqs <- function(class=NULL, complement=FALSE, not_computed=FALSE, after_id=0)
{
	query = paste("select seq_id as seq from sequence_class");
	
	query = paste(query, " where id >", after_id);
	
	#kdrew: if class is specified
	if(!is.null(class))
	{
		#kdrew: if complement flag is not set just get the records in class
		if(FALSE == complement)
		{
			query = paste(query, " and class =", class);
		}
		#kdrew: if complement flag is set, get records not in class
		else
		{
			query = paste(query, " and class <>", class);
		}
		#kdrew: only get sequences which have not been completed
		if(TRUE == not_computed)
		{	
			query = paste(query, " and prob_computed = 0 ");
		}
	}
	sequence_class = hpf_qry(query);
	return(sequence_class);
}

#kdrew: return the classes sequences are divided into  
getClasses <- function(class=NULL, complement=FALSE)
{
	query = paste("select DISTINCT class from sequence_class");
	if(!is.null(class))
	{
		#kdrew: if complement flag is not set just get the records in class
		if(FALSE == complement)
		{
			query = paste(query, " where class =", class);
		}
		#kdrew: if complement flag is set, get records not in class
		else
		{
			query = paste(query, " where class <>", class);
		}
	}
	classes = hpf_qry(query);
	return(classes[,1]);
}

is_mutually_exclusive <- function(go_terms)
{
	query = paste("select * from graph_path ");
	query = paste(query, " where ");
	i = 0;
	
	lapply(go_terms, conditional_me <- function(gt1)
			{
				lapply(go_terms, conditional_me2 <- function(gt2)
						{
							if(gt1 != gt2)
							{
								if(0 < i)
								{
									query <<- paste(query, " or ");
								}
								query <<- paste(query, "term1_id = ",gt1,"and term2_id = ",gt2);
								i <<- i + 1;
							}
						});
			});
	
	print(query);
	
	mygo_qry(query);
	
}


