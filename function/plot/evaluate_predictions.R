
source("sql_setup.R");
source("go_hpf_utils.R");
source("create_dbs.R");
source("compute_function_probabilities.R");
source("compute_function_bayes.R");


#create temporary table mf_baseline2 (select *, (select count(*) + 1 from mf_baseline as mb2 where mb2.prob > mb.prob) as rank from mf_baseline as mb);
#create temporary table 40761_eval (select *, (base_rank - rank) as rank_diff, (prob - base_prob) as prob_diff from (select s2mp.*, (select count(*) +1 from seq2mf_prob_tmp where seq_id = s2mp.seq_id and is_max = 1 and prob > s2mp.prob) as rank, mb2.term_count, mb2.prob as base_prob, mb2.rank as base_rank from mf_baseline as mb2, seq2mf_prob_tmp as s2mp where s2mp.mf_id = mb2.term_id and s2mp.seq_id = 40761 and is_max = 1 order by prob DESC) as tmp);

evaluate_benchmark <- function(benchmark_db, sf_data_table="mcm_data", rank_thres=5, scale_quantity=.8, scale_multiple=FALSE, training_set_suffix, evaluation_table, bayes=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA','NR'))
{
	if(bayes)
	{
		query = paste("create table if not exists ",evaluation_table," like hpf.evaluation_results_bayes");
		hpf_qry(query);
	}
	else
	{
		query = paste("create table if not exists ",evaluation_table," like hpf.evaluation_results");
		hpf_qry(query);
	}

	query = paste("select distinct domain_sequence_key from ",benchmark_db,".",sf_data_table," where domain_sequence_key not in (select distinct seq_id from ",evaluation_table,")", sep="");

	seqs = hpf_qry(query);

	apply(seqs,1,function(seq)
	{
		#kdrew: make predictions and store all mfs
		evaluate_predictions(seq, testSet=benchmark_db, sum_sf=TRUE, evaluation_table=evaluation_table, allmfs=TRUE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, scale_multiple=scale_multiple, training_set_suffix=training_set_suffix,bayes=bayes, code=code);
	});

}



#kdrew: takes in sequence, computes probabilities, ranks and then compares predictor sets
#kdrew: make predictions, default to store only known/true mfs , ie if allmfs is true then all results are stored
evaluate_predictions <- function(seq_id, testSet="goSet", sum_sf=FALSE, evaluation_table="evaluation_results", allmfs=FALSE, sf_data_table,  rank_thres, scale_quantity, scale_multiple, training_set_suffix, bayes=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA','NR'))
{
	class=0;
	class_comp = getClasses(class,complement=TRUE);

	create_seq2mf_prob("seq2mf_prob_tmp_pls", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_pl", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_ps", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_ls", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_p", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_l", drop=TRUE);
	create_seq2mf_prob("seq2mf_prob_tmp_s", drop=TRUE);

	if(bayes)
	{
		store_table = "bayes_log_ratio_tmp";
		filtered_table = "bayes_log_ratio_filtered";

		bayes_classifer_continous_sf_sum2(seq_id=seq_id, process=TRUE, localization=TRUE, structure=TRUE, 
						store_table=store_table, testSet=testSet, sf_data_table=sf_data_table, 
						rank_thres=rank_thres, scale_quantity=scale_quantity, 
						scale_multiple=scale_multiple, training_set_suffix=training_set_suffix, code=code);

		filter_bayes_results(filtered_table=filtered_table, results_table=store_table, training_set_suffix=training_set_suffix);

		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_pls", BP=TRUE, CC=TRUE, SF=TRUE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_pls")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_pl", BP=TRUE, CC=TRUE, SF=FALSE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_pl")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_ps", BP=TRUE, CC=FALSE, SF=TRUE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_ps")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_ls", BP=FALSE, CC=TRUE, SF=TRUE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_ls")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_p", BP=TRUE, CC=FALSE, SF=FALSE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_p")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_l", BP=FALSE, CC=TRUE, SF=FALSE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_l")
		store_bayes_results(seq_id=seq_id, store_table=filtered_table, table="seq2mf_prob_tmp_s", BP=FALSE, CC=FALSE, SF=TRUE);
		include_missed_mf_bayes(seq=seq_id, table="seq2mf_prob_tmp_s")

		true = compile_results_bayes(seq_id, evaluation_table=evaluation_table);
	}
	else
	{
		#kdrew: create and populate seq2mf_prob tables
		computeProbabilities(seq_id,table="seq2mf_prob_tmp_pls", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,structure=FALSE,table="seq2mf_prob_tmp_pl", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,localization=FALSE,table="seq2mf_prob_tmp_ps", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,process=FALSE,table="seq2mf_prob_tmp_ls", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,localization=FALSE,structure=FALSE,table="seq2mf_prob_tmp_p", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,process=FALSE, structure=FALSE,table="seq2mf_prob_tmp_l", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);

		computeProbabilities(seq_id,process=FALSE, localization=FALSE,table="seq2mf_prob_tmp_s", testSet=testSet, sum_sf=sum_sf, combination=FALSE, sf_data_table=sf_data_table, rank_thres=rank_thres, scale_quantity=scale_quantity, training_set_suffix=training_set_suffix,bayes=bayes);


		#kdrew: determine ranks for functions under proc, loc and sf
		create_rank_table(seq_id,"pls")
		#kdrew: determine ranks for functions under proc and loc
		create_rank_table(seq_id,"pl")
		#kdrew: determine ranks for functions under proc and sf
		create_rank_table(seq_id,"ps")
		#kdrew: determine ranks for functions under loc and sf
		create_rank_table(seq_id,"ls")
		#kdrew: determine ranks for functions under proc
		create_rank_table(seq_id,"p")
		#kdrew: determine ranks for functions under loc
		create_rank_table(seq_id,"l")
		#kdrew: determine ranks for functions under sf
		create_rank_table(seq_id,"s")

		true = compile_results(seq_id, evaluation_table=evaluation_table);

		drop_rank_table(seq_id,"s")
		drop_rank_table(seq_id,"l")
		drop_rank_table(seq_id,"p")
		drop_rank_table(seq_id,"pl")
		drop_rank_table(seq_id,"ls")
		drop_rank_table(seq_id,"ps")
		drop_rank_table(seq_id,"pls")
	}
	
	return(true);
}

#select * from 40761_eval as ev LEFT JOIN mfs on (ev.mf_id = mfs.mf_id) where mfs.mf_id is not null limit 15;
compile_results <- function(seq_id, db_suffix = "", evaluation_table, bayes=FALSE)
{
	#query = paste("select ev.id, ev.seq_id, ev.mf_id, ev.bp_id, ev.cc_id, ev.sf_id, ev.prob, ev.rank, ev.base_prob, ev.base_rank, ev.rank_diff");
	#query = paste(query, " from ",seq_id,"_eval_",db_suffix," as ev LEFT JOIN mfs on (ev.mf_id = mfs.mf_id) where mfs.mf_id is not null", sep="");

	ev_p = paste(seq_id,"_eval_p", sep="");
	ev_l = paste(seq_id,"_eval_l", sep="");
	ev_s = paste(seq_id,"_eval_s", sep="");
	ev_pl = paste(seq_id,"_eval_pl", sep="");
	ev_ps = paste(seq_id,"_eval_ps", sep="");
	ev_ls = paste(seq_id,"_eval_ls", sep="");
	ev_pls = paste(seq_id,"_eval_pls", sep="");

	query = paste("insert into ", evaluation_table);

	query = paste(query, " select ev_pls.seq_id, ev_pls.mf_id, t.name ");

	query = paste(query, ", ev_p.prob as p_prob");
	query = paste(query, ", ev_l.prob as l_prob");
	query = paste(query, ", ev_s.prob as s_prob");
	query = paste(query, ", ev_pl.prob as pl_prob");
	query = paste(query, ", ev_ps.prob as ps_prob");
	query = paste(query, ", ev_ls.prob as ls_prob");
	query = paste(query, ", ev_pls.prob as pls_prob");
	query = paste(query, ", ev_pls.base_prob");

	query = paste(query, ", ev_p.rank as p_rank");
	query = paste(query, ", ev_l.rank as l_rank");
	query = paste(query, ", ev_s.rank as s_rank");
	query = paste(query, ", ev_pl.rank as pl_rank");
	query = paste(query, ", ev_ps.rank as ps_rank");
	query = paste(query, ", ev_ls.rank as ls_rank");
	query = paste(query, ", ev_pls.rank as pls_rank");
	query = paste(query, ", ev_pls.base_rank ");

	query = paste(query, " from ");
	query = paste(query, ev_pls, "as ev_pls");

	query = paste(query, ",", ev_p, " as ev_p ");
	query = paste(query, ",", ev_l, " as ev_l ");
	query = paste(query, ",", ev_s, " as ev_s ");
	query = paste(query, ",", ev_pl," as ev_pl");
	query = paste(query, ",", ev_ps," as ev_ps");
	query = paste(query, ",", ev_ls," as ev_ls");

	query = paste(query, ", mygo.term as t");

	query = paste(query, " where ev_pls.mf_id is not null");

	query = paste(query, " and ev_p.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_l.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_s.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_pl.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_ps.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_ls.mf_id = ev_pls.mf_id");

	query = paste(query, " and t.id = ev_pls.mf_id");

	query = paste(query, " order by ev_pls.base_rank ");

	print(query);
	hpf_qry(query);
}

compile_results_bayes <- function(seq_id, evaluation_table=evaluation_table)
{
	ev_p = paste("seq2mf_prob_tmp_p");
	ev_l = paste("seq2mf_prob_tmp_l");
	ev_s = paste("seq2mf_prob_tmp_s");
	ev_pl = paste("seq2mf_prob_tmp_pl");
	ev_ps = paste("seq2mf_prob_tmp_ps");
	ev_ls = paste("seq2mf_prob_tmp_ls");
	ev_pls = paste("seq2mf_prob_tmp_pls");

	query = paste("insert into ", evaluation_table);

	query = paste(query, " select ev_pls.seq_id, ev_pls.mf_id, t.name ");

	query = paste(query, ", ev_p.log_likelihood_ratio as p_llr");
	query = paste(query, ", ev_l.log_likelihood_ratio as l_llr");
	query = paste(query, ", ev_s.log_likelihood_ratio as s_llr");
	query = paste(query, ", ev_pl.log_likelihood_ratio as pl_llr");
	query = paste(query, ", ev_ps.log_likelihood_ratio as ps_llr");
	query = paste(query, ", ev_ls.log_likelihood_ratio as ls_llr");
	query = paste(query, ", ev_pls.log_likelihood_ratio as pls_llr");
	query = paste(query, ", log((mfb.prob+",TINY_NUM,")/(1-mfb.prob+",TINY_NUM,")) as base_llr");

	query = paste(query, " from ");
	query = paste(query, ev_pls, "as ev_pls");
	query = paste(query, ",", ev_p, " as ev_p ");
	query = paste(query, ",", ev_l, " as ev_l ");
	query = paste(query, ",", ev_s, " as ev_s ");
	query = paste(query, ",", ev_pl," as ev_pl");
	query = paste(query, ",", ev_ps," as ev_ps");
	query = paste(query, ",", ev_ls," as ev_ls");

	query = paste(query, ", mygo.term as t");
	query = paste(query, ", mf_baseline as mfb");

	query = paste(query, " where ev_pls.mf_id is not null");

	query = paste(query, " and ev_p.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_l.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_s.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_pl.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_ps.mf_id = ev_pls.mf_id");
	query = paste(query, " and ev_ls.mf_id = ev_pls.mf_id");

	query = paste(query, " and t.id = ev_pls.mf_id");
	query = paste(query, " and t.id = mfb.mf_id");

	print(query);
	hpf_qry(query);
}

drop_rank_table <- function(seq_id, db_suffix = "")
{
	drop_query = paste("drop table if exists ", seq_id,"_eval_",db_suffix,sep="");
	print(drop_query);
	hpf_qry(drop_query);
}

create_rank_table <- function(seq_id, db_suffix = "")
{
#select @rank:=0; create temporary table 40761_eval_ls (select *, (base_rank - rank) as rank_diff, (prob - base_prob) as prob_diff from (select s2mp.*, @rank:=@rank+1 as rank, mb2.term_count, mb2.prob as base_prob, mb2.rank as base_rank from mf_baseline as mb2, seq2mf_prob_tmp_ls as s2mp where s2mp.mf_id = mb2.term_id and s2mp.seq_id = 40761 and is_max = 1 order by prob DESC) as tmp);

	#drop_query = paste("drop table if exists ", seq_id,"_eval_",db_suffix,sep="");

	rank_query = "select @rank:=0";

	query = paste("create temporary table ", seq_id,"_eval_",db_suffix, sep=""); 
	#kdrew: inserted rank_diff_rank and prob_diff_rank
	query = paste(query, " (select *, (base_rank - rank) as rank_diff, (prob - base_prob) as prob_diff, 0 as rank_diff_rank, 0 as prob_diff_rank");
	query = paste(query, " from (select s2mp.*, @rank:=@rank+1 as rank, mb2.term_count, mb2.prob as base_prob, mb2.rank as base_rank ");
	query = paste(query, " from mf_baseline as mb2, seq2mf_prob_tmp_",db_suffix,sep=''); 
	query = paste(query, " as s2mp where s2mp.mf_id = mb2.term_id and s2mp.seq_id = ", seq_id);
	query = paste(query, " and is_max = 1 ");
	query = paste(query, " order by prob DESC) as tmp)");

	alter_query = paste("alter table ", seq_id,"_eval_",db_suffix, sep="");
	alter_query = paste(alter_query, " add index(mf_id)");

	#print(drop_query);
	print(rank_query);
	print(query);
	print(alter_query);

	#hpf_qry(drop_query);
	drop_rank_table(seq_id,db_suffix);

	hpf_qry(rank_query);
	hpf_qry(query);
	hpf_qry(alter_query);

	drop_rank_query = paste("drop table if exists rank_diff_tmp");
	drop_prob_query = paste("drop table if exists prob_diff_tmp");

	#kdrew: TODO parameterize rank <= 100
	rank_diff_query = paste("create temporary table rank_diff_tmp select mf_id, @rank := @rank + 1 as rank_diff_rank from ", seq_id,"_eval_",db_suffix," where rank <= 100 order by rank_diff DESC", sep="");
	prob_diff_query = paste("create temporary table prob_diff_tmp select mf_id, @rank := @rank + 1 as prob_diff_rank from ", seq_id,"_eval_",db_suffix," where rank <= 100 order by prob_diff DESC", sep="");

	print(drop_rank_query);
	print(drop_prob_query);
	hpf_qry(drop_rank_query);
	hpf_qry(drop_prob_query);

	print(rank_query);
	print(rank_diff_query);

	hpf_qry(rank_query);
	hpf_qry(rank_diff_query);

	print(rank_query);
	print(prob_diff_query);

	hpf_qry(rank_query);
	hpf_qry(prob_diff_query);
	
	update_query = paste("update ",seq_id,"_eval_",db_suffix," as ev, rank_diff_tmp as rdt set ev.rank_diff_rank = rdt.rank_diff_rank where ev.mf_id = rdt.mf_id", sep="");
	update2_query = paste("update ",seq_id,"_eval_",db_suffix," as ev, prob_diff_tmp as pdt set ev.prob_diff_rank = pdt.prob_diff_rank where ev.mf_id = pdt.mf_id", sep="");

	print(update_query);
	print(update2_query);

	hpf_qry(update_query);
	hpf_qry(update2_query);
}


##drop table if exists mfs; create temporary table mfs (select DISTINCT t.id as mf_id  from seq2go as s2g, sequence_class as sc, mygo.gene_product_seq as gps, mygo.term as t  where s2g.gene_product_id = gps.gene_product_id and gps.seq_id = sc.seq_id   and t.id = s2g.term1_id and t.term_type = 'molecular_function'   and s2g.code IN ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA')  and sc.seq_id = 40761);
#
##kdrew: this creates a table which lists ids for all molecular functions in the GO database (not just the ones annotated to the testSet)
#create_mfs <- function(seq_id, allmfs = FALSE)
#{
#		query_drop = paste("drop table if exists mfs");
#
#		query = paste("create temporary table mfs (");
#		query = paste(query, "select DISTINCT t.id as mf_id ");
#		query = paste(query, " from mygo.term as t ");
#		query = paste(query, " where t.term_type = 'molecular_function'  )"); 
#
#		print(query_drop);
#		print(query);
#
#		hpf_qry(query_drop);
#		hpf_qry(query);
#}


randomize_seqs <- function(seq_table, rand_seq_table="rand_seq_map", offset=500)
{
	
	drop_query = paste("drop table if exists", rand_seq_table);
	create_query = paste("create table ",rand_seq_table, " select seq_id, -1 as id_rand, -1 as rand_seq_id from ",seq_table);
	key_query = paste("alter table ",rand_seq_table," ADD id INT AUTO_INCREMENT PRIMARY KEY");

	#kdrew: old version of randomizing sequences, didn't care about duplicate sequences
	#rand_query = paste("update ",rand_seq_table, "set rand_seq_id = (select seq_id from ",seq_table, " order by rand(",offset,") limit 1)");

	#select @rand_cnt := 0;
	var_query = paste("select @rand_cnt := 0");
	#update seqs set id_rand = (@rand_cnt := @rand_cnt + 1) order by rand(10);
	rand_query = paste("update ", rand_seq_table, " set id_rand = (@rand_cnt := @rand_cnt + 1) order by rand(",offset,")");
	#update seqs as s1, seqs as s2 set s1.seq_rand = s2.seq_id where s1.id_rand = s2.id;
	update_query = paste("update ", rand_seq_table, " as s1, ", rand_seq_table, " as s2 set s1.rand_seq_id = s2.seq_id where s1.id_rand = s2.id");

	print(drop_query);
	print(create_query);
	print(key_query);
	print(var_query);
	print(rand_query);
	print(update_query);

	hpf_qry(drop_query);
	hpf_qry(create_query);
	hpf_qry(key_query);
	hpf_qry(var_query);
	hpf_qry(rand_query);
	hpf_qry(update_query);

	alter_query = paste("alter table ",rand_seq_table," add index(seq_id)");
	hpf_qry(alter_query);
	alter_query = paste("alter table ",rand_seq_table," add index(rand_seq_id)");
	hpf_qry(alter_query);
}

############
#kdrew: takes in a list of sequences in a table and computes the number of correct functions and the number of all functions for given threshold values
#eval_table is table with computed probabilities for all functions for given sequences, functionFrom is a table of true functions matched to sequences, predictors is the type of evidence to be evaluated
#random if true will shuffle the function lables, offset is a seed for shuffling lables, score_table is table to deposit results
############
correctOverAllPredictions <- function(seq_table, prob_score_thres_list = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9),rank_score_thres_list = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9), rank_thres = 100, eval_table = "evaluation_results_PDB", functionFrom = "scopfold", predictors="s", random=FALSE, offset=500, score_table, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA', 'IC','IEP','NAS','ND','NR','RCA'))
{

	rand_seq_table="rand_seq_map";

	drop_query = paste("drop table if exists evalPred");
	create_query = paste("create temporary table evalPred select ery.seq_id, mf_id, ",predictors,"_prob, base_prob, ",predictors,"_rank, base_rank", sep="");
	create_query = paste(create_query, " , if(",predictors,"_prob = base_prob, 0, (",predictors,"_prob - base_prob)/(1.0-base_prob)) as prob_score, (base_rank - ",predictors,"_rank)/base_rank as rank_score ", sep="");
	create_query = paste(create_query, " from ",eval_table," as ery, ",seq_table," as sit where ery.seq_id = sit.seq_id");

	print(drop_query);
	hpf_qry(drop_query);
	print(create_query);
	hpf_qry(create_query);

	alter_query = paste("alter table evalPred add index(seq_id)");
	hpf_qry(alter_query);
	alter_query = paste("alter table evalPred add index(mf_id)");
	hpf_qry(alter_query);

	getFunction(code = code, functionFrom = functionFrom, create_table=TRUE, seq_table=seq_table);

	for(i in seq(length(prob_score_thres_list)))
	{
		for(j in seq(length(rank_score_thres_list)))
		{
			is_random = 0;
			if(random)
			{
				randomize_seqs(seq_table,rand_seq_table=rand_seq_table, offset=offset);

				count_notNULL_query = paste("select count(*) from evalPred as ep left outer join mfs on (ep.mf_id = mfs.mf_id and (select rand_seq_id from ",rand_seq_table," where seq_id = ep.seq_id) = mfs.seq_id) ");
				count_notNULL_query = paste(count_notNULL_query, " where ep.prob_score > ",prob_score_thres_list[i]," and ep.rank_score > ",rank_score_thres_list[j]," and mfs.seq_id is not null");

				count_query = paste("select count(*) from evalPred as ep left outer join mfs on (ep.mf_id = mfs.mf_id and (select rand_seq_id from ",rand_seq_table," where seq_id = ep.seq_id) = mfs.seq_id) ");
				count_query = paste(count_query, " where ep.prob_score > ",prob_score_thres_list[i]," and ep.rank_score > ",rank_score_thres_list[j]);

				is_random = 1;
				
			}
			else
			{

				count_notNULL_query = paste("select count(*) from evalPred as ep left outer join mfs on (ep.mf_id = mfs.mf_id and ep.seq_id = mfs.seq_id) ");
				count_notNULL_query = paste(count_notNULL_query, " where ep.prob_score > ",prob_score_thres_list[i]," and ep.rank_score > ",rank_score_thres_list[j]," and mfs.seq_id is not null");

				count_query = paste("select count(*) from evalPred as ep left outer join mfs on (ep.mf_id = mfs.mf_id and ep.seq_id = mfs.seq_id) ");
				count_query = paste(count_query, " where ep.prob_score > ",prob_score_thres_list[i]," and ep.rank_score > ",rank_score_thres_list[j]);
			}

			print(count_query);
			full_count = hpf_qry(count_query);
			print(count_notNULL_query);
			real_count = hpf_qry(count_notNULL_query);

			bookeeping_query = paste("insert into ",score_table,"( eval_table, functionFrom, predictors, prob_score_thres, rank_score_thres, rank_thres, real_count, full_count, random, offset) VALUES (\"",eval_table,"\",\"",functionFrom,"\",\"",predictors,"\",",prob_score_thres_list[i],",",rank_score_thres_list[j],",",rank_thres,",",real_count,",",full_count,",\"",is_random,"\",",offset,")");
			print(bookeeping_query);
			hpf_qry(bookeeping_query);
		}
	}

}

