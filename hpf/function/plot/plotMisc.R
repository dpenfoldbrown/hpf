
source("sql_setup.R");
source("go_hpf_utils.R");
#source("evaluate_predictions.R");

functionCount <- function(seq_table, eval_table, functionFrom, predictors="s", break_points = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.0), code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA', 'IC','IEP','NAS','ND','NR','RCA'))
{
	func_cnt_vector = c();
	seq_cnt_vector = c();
	func_cnt_nrm_vector = c();

	getFunction(code=code, functionFrom = functionFrom, create_table=TRUE, seq_table=seq_table);

	drop_query = paste("drop table if exists correctFunctions"); 
	print(drop_query);
	hpf_qry(drop_query)

	query = paste("create temporary table correctFunctions select ery.domain_sequence_key, ery.mf_acc, ery.name, ery.p_prob, ery.s_prob, ery.ps_prob, ery.base_prob, mfs.mf_acc as mfs_mf_acc from ", eval_table, " as ery left outer join mfs on (ery.domain_sequence_key = mfs.seq_id and ery.mf_acc = mfs.mf_acc)");
	print(query);
	hpf_qry(query)

	alter_query = paste("alter table correctFunctions add index(seq_id)");
	print(alter_query);
	hpf_qry(alter_query)
	alter_query = paste("alter table correctFunctions add index(mf_acc)");
	print(alter_query);
	hpf_qry(alter_query)

	for(i in seq(length(break_points)))
	{

		query = paste("select count(distinct mf_acc, seq_id) from correctFunctions where mfs_mf_acc is not NULL and base_prob > ", break_points[i]);
		print(query);
		function_count = hpf_qry(query)
		print(paste("function_count:",function_count));

		query = paste("select count(distinct seq_id) from correctFunctions where mfs_mf_acc is not NULL and base_prob > ", break_points[i]);
		print(query);
		seq_count = hpf_qry(query)
		print(paste("seq_count:",seq_count));

		func_cnt_vector <- c(func_cnt_vector, function_count[1,1]);
		seq_cnt_vector <- c(seq_cnt_vector, seq_count[1,1]);
		func_cnt_nrm_vector <- c(func_cnt_nrm_vector, function_count[1,1]/seq_count[1,1]);
	}

		print(func_cnt_vector);
		print(seq_cnt_vector);
		print(func_cnt_nrm_vector);

		return(cbind(func_cnt_vector, seq_cnt_vector, func_cnt_nrm_vector));
}

since_solved_analysis <- function(seq_table="", benchmark="scop_benchmark", sf_data_table="mcm_data", since_solved_table="hddb.fold_sccs_20071107", color1="grey", color2="yellow", plot=TRUE)
{
	ret = list();

	#scop_mcm_scores_true = hpf_qry("select probability from hddb.fold_sccs_20071107 as fs, scop_benchmark.mcm_data as md where md.parent_sequence_key = fs.sequence_key and md.sccs = fs.sccs")

	all_query = paste("select probability from ", since_solved_table, " as fs, ");
	all_query = paste(all_query, benchmark,".", sf_data_table, " as md ", sep="");
	if(seq_table!="")
	{ all_query = paste(all_query, ",", seq_table, " as seq_table ", sep=""); }
	all_query = paste(all_query, " where md.parent_sequence_key = fs.sequence_key");
	if(seq_table!="")
	{ all_query = paste(all_query, " and seq_table.seq_id = md.domain_sequence_key "); }
	print(all_query);
	mcm_scores = hpf_qry(all_query);

	true_query = paste("select probability from ", since_solved_table, " as fs, ");
	true_query = paste(true_query, benchmark,".", sf_data_table, " as md ", sep="");

	if(seq_table!="")
	{ true_query = paste(true_query, ",",  seq_table, " as seq_table ", sep=""); }
	true_query = paste(true_query, " where md.parent_sequence_key = fs.sequence_key and md.sccs=fs.sccs");

	if(seq_table!="")
	{ true_query = paste(true_query, " and seq_table.seq_id = md.domain_sequence_key "); }

	print(true_query);
	mcm_scores_true = hpf_qry(true_query);

	incorrect = mcm_scores[,1][!mcm_scores[,1] %in% mcm_scores_true[,1]]; 
	ks_true_v_incorrect = ks.test(mcm_scores_true[,1], incorrect)
	pval = ks_true_v_incorrect$p.value
	dstat = ks_true_v_incorrect$statistic

	if(plot)
	{
		since_solved_plotHist(benchmark=benchmark, scores=mcm_scores, scores_true=mcm_scores_true, pval=pval, dstat=dstat, color1=color1, color2=color2);
	}
	
	ret$mcm_scores = mcm_scores[,1];
	ret$mcm_scores_true = mcm_scores_true[,1];

	return(ret);
}

since_solved_plotHist <-function(benchmark, scores, scores_true, pval, dstat, color1="grey", color2="yellow")
{
	pval = format(pval,scientific=TRUE)
	dstat = format(dstat,digits=2)
	quartz();
	hist_tmp = hist(as.matrix(scores), xlim=c(0,1), xlab="mcm_scores", ylab="count", col=color1, main=paste("Histogram of",benchmark,":",length(scores_true[,1])," true /",length(scores[,1])," total"), sub=paste(color1,"= all mcm scores, ",color2,"= correct based on since solved, KS-test: D=",dstat,"p-value=", pval), breaks=seq(0,1,by=.05), cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2)
	hist_true_tmp = hist(as.matrix(scores_true), breaks=hist_tmp$breaks, col=color2, add=TRUE)
	labels=hist_true_tmp$counts/hist_tmp$counts;
	label_nans = is.nan(labels)
	labels=format(labels,digits=1)
	labels[label_nans] = ""
	#labels[labels%in%"NaN"] = ""
	text(hist_tmp$mids,hist_tmp$count,labels=labels, pos=3, cex=2)
}

ks_test_analysis <- function(x,y=x)
{
	ks_ret = list();
	x_incorrect = x$mcm_scores[!x$mcm_scores %in% x$mcm_scores_true]; 
	y_incorrect = y$mcm_scores[!y$mcm_scores %in% y$mcm_scores_true]; 

	ks_ret$all_v_all = ks.test(x$mcm_scores,y$mcm_scores)
	ks_ret$true_v_true = ks.test(x$mcm_scores_true,y$mcm_scores_true)
	ks_ret$all_v_true = ks.test(x$mcm_scores,y$mcm_scores_true)
	ks_ret$true_v_all = ks.test(x$mcm_scores_true,y$mcm_scores)

	ks_ret$incorrect_v_incorrect = ks.test(x_incorrect, y_incorrect);
	ks_ret$incorrect_v_true = ks.test(x_incorrect, y$mcm_scores_true);
	ks_ret$true_v_incorrect = ks.test(x$mcm_scores_true,y_incorrect)

	return(ks_ret);
}

length_plotHist <- function(data_table_name, field_name="domain_sequence_key", color="grey")
{

	len_qry = paste("select dl.len from hddb.ddb_length as dl where dl.sequence_key in (");
	len_qry = paste(len_qry, " select distinct ", field_name, " from ", data_table_name);
	len_qry = paste(len_qry, " )");

	print(len_qry);
	lengths = hpf_qry(len_qry);

	quartz();
	hist_tmp = hist(as.matrix(lengths), xlab="lengths", ylab="count", col=color, main=paste("histogram of lengths",data_table_name), breaks=c(seq(20,2000,by=10)));
	#hist_tmp = hist(as.matrix(lengths), xlab="lengths", ylab="count", col=color, main=paste("histogram of lengths",data_table_name));

}

#kdrew: for a given benchmark and sequences in the benchmark create a histogram of the number of annotations to a protein 
functionResultsCount_plotHist <- function(seq_table, eval_table, predictors="s", llr_min, llr_max, base_max, color="grey", break_points="Sturges")
{
	#select base_llr  from swissprot_mf_seqs_sinceSolved as ss, hpf_results.swissprot_bayes_results as sr where sr.sequence_key = ss.seq_id and s_llr <= 7 and s_llr > 2 and base_llr < 0 limit 5

	query = paste("select base_llr ");
	query = paste(query, " from ", seq_table, " as st, ", eval_table," as et ");
	query = paste(query, " where st.seq_id = et.sequence_key ");
	query = paste(query, " and ",predictors,"_llr <= ", llr_max," and ",predictors,"_llr > ",llr_min," and base_llr < ", base_max, sep="");

	print(query);
	base_llrs = as.matrix(hpf_qry(query));

	quartz();
	hist_tmp = hist(base_llrs, xlab="base log likelihood ratio", ylab="frequency", col=color, main=paste("histogram of function base llr :",seq_table),breaks=break_points);
}

#kdrew: for a given benchmark and sequences in the benchmark create a histogram of the number of annotations to a protein 
functionCount_plotHist <- function(seq_table, data_table="mcm_data", benchmark_db, color="grey", breaks=c(seq(0,20,by=1)))
{
	#select bs.seq_id, count(distinct t.id) as cnt from bacteria_sp_mf_seqs as bs, mcm_data as md, function as f, mygo.term as t where bs.seq_id = md.domain_sequence_key and md.parent_sequence_key = f.sequence_key and f.acc = t.acc and t.term_type = "molecular_function" group by bs.seq_id
	
	cnt_query = paste("select count(distinct t.id) as cnt ");
	cnt_query = paste(cnt_query, " from ", seq_table, " as st, ", benchmark_db,".",data_table, " as dt, ",benchmark_db,".function as f, mygo.term as t", sep=""); 
	cnt_query = paste(cnt_query, " where st.seq_id = dt.domain_sequence_key and dt.parent_sequence_key = f.sequence_key and f.acc = t.acc and t.term_type = 'molecular_function'");
	cnt_query = paste(cnt_query, "  group by st.seq_id");

	print(cnt_query);
	counts = as.matrix(hpf_qry(cnt_query));

	quartz();
	hist_tmp = hist(counts, xlab="function counts", ylab="frequency", col=color, main=paste("histogram of function counts ",benchmark_db, ":",seq_table),breaks=breaks);
	#hist_tmp = hist(as.matrix(lengths), xlab="lengths", ylab="count", col=color, main=paste("histogram of lengths",data_table_name));
}

#funcPred_plotHist(seq_table="swissprot_benchmark.swissprot_mf_seqs_1000", eval_table="hpf_results.swissprot_bayes_results_mf_1000", functionFrom="swissprot_benchmark")
funcPred_plotHist <- function(seq_table, eval_table, functionFrom, predictors="s", break_points="Sturges", ylim=NULL, predictor_min = 0, base_max = 0, base_min = -99, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA', 'IC','IEP','NAS','ND','NR','RCA'), bayes=TRUE, color1="grey", color2="seagreen", title="", percent_offset=.03, yield_points=c(3,8), random=FALSE, rand_seq_table="rand_seq_map", offset=500)
{
	ret_val = list();
	getFunction(code=code, functionFrom = functionFrom, create_table=TRUE, seq_table=seq_table);

	#swissprot_correct = hpf_qry("select ery.domain_sequence_key, ery.mf_acc, ery.s_llr, ery.base_llr, mfs.mf_acc as mfs_mf_acc   from  evaluation_results_swissprot_bayes8_sq9_smulti2_part  as ery, mfs, swissprot_benchmark.swissprot_mf_seqs as seq_table  where ery.domain_sequence_key = mfs.seq_id and ery.mf_acc = mfs.mf_acc and seq_table.seq_id = ery.domain_sequence_key   and ery.base_llr <=  0  and ery.base_llr >=  -99 and ery.s_llr >= 0")
	#swissprot_all = hpf_qry("select ery.domain_sequence_key, ery.mf_acc, ery.s_llr, ery.base_llr from  evaluation_results_swissprot_bayes8_sq9_smulti2_part  as ery, swissprot_benchmark.swissprot_mf_seqs as seq_table  where seq_table.seq_id = ery.domain_sequence_key   and ery.base_llr <=  0  and ery.base_llr >=  -99 and ery.s_llr >= 0")

		if(bayes)
		{
			#cor_query = paste(" select distinct ery.domain_sequence_key, ery.mf_acc, ery.",predictors,"_llr as llr", sep="");
			#all_query = paste(" select distinct ery.domain_sequence_key, ery.mf_acc, ery.",predictors,"_llr as llr", sep="");
			cor_query = paste(" select ery.",predictors,"_llr as llr", sep="");
			all_query = paste(" select ery.",predictors,"_llr as llr", sep="");
		}
		else
		{
			#cor_query = paste(" select ery.domain_sequence_key, ery.mf_acc, ery.",predictors,"_prob as prob", sep="");
			#all_query = paste(" select ery.domain_sequence_key, ery.mf_acc, ery.",predictors,"_prob as prob", sep="");
			cor_query = paste(" select ery.",predictors,"_prob as prob", sep="");
			all_query = paste(" select ery.",predictors,"_prob as prob", sep="");
		}

		cor_query = paste(cor_query, " from ",eval_table, " as ery ");
		all_query = paste(all_query, " from ",eval_table, " as ery ");

		#cor_query = paste(cor_query, " , ", seq_table," as seq_table");
		#all_query = paste(all_query, " , ", seq_table," as seq_table");
		cor_query = paste(cor_query, " , mfs ");
		cor_query = paste(cor_query, " , ", seq_table," as seq_table, mygo.term as t");
		all_query = paste(all_query, " , ", seq_table," as seq_table");

		if(random)
		{
			randomize_seqs(seq_table,rand_seq_table=rand_seq_table, offset=offset);

			cor_query = paste(cor_query, ", ", rand_seq_table," as rand_seq_table ");
			cor_query = paste(cor_query," where seq_table.seq_id = rand_seq_table.seq_id and seq_table.seq_id = ery.domain_sequence_key and rand_seq_table.rand_seq_id = mfs.seq_id and ery.mf_acc = t.acc and t.id = mfs.mf_acc");
		}
		else
		{
			#cor_query = paste(cor_query, " where ery.domain_sequence_key = mfs.seq_id and ery.mf_acc = mfs.mf_acc and seq_table.seq_id = ery.domain_sequence_key ");
			#all_query = paste(all_query, " where seq_table.seq_id = ery.domain_sequence_key ");
			cor_query = paste(cor_query, " where ery.domain_sequence_key = mfs.seq_id and ery.mf_acc = t.acc and t.id = mfs.mf_acc and seq_table.seq_id = ery.domain_sequence_key");
		}

		all_query = paste(all_query, " where seq_table.seq_id = ery.domain_sequence_key");

		if(bayes)
		{
			cor_query = paste(cor_query, " and ery.base_llr <= ",base_max, " and ery.base_llr >= ",base_min);
			cor_query = paste(cor_query, " and ery.",predictors,"_llr >= ", predictor_min, sep="");

			all_query = paste(all_query, " and ery.base_llr <= ",base_max, " and ery.base_llr >= ",base_min);
			all_query = paste(all_query, " and ery.",predictors,"_llr >= ", predictor_min, sep="");
		}
		else
		{
			cor_query = paste(cor_query, " and ery.base_prob <= ",base_max, " and ery.base_prob >= ",base_min);
			cor_query = paste(cor_query, " and ery.",predictors,"_prob >= ", predictor_min, sep="");

			all_query = paste(all_query, " and ery.base_prob <= ",base_max, " and ery.base_prob >= ",base_min);
			all_query = paste(all_query, " and ery.",predictors,"_prob >= ", predictor_min, sep="");
		}

		print(cor_query);
		correct = hpf_qry(cor_query)

		print(all_query);
		all = hpf_qry(all_query)

		print("done with all_query")

		if(title == "")
		{
			title = paste("Histogram of Function Prediction for ", functionFrom, ": ",predictors, " predictors"); 
		}
		print("after with title")

		domain_cnt = getDomainCount(eval_table,seq_table);

		sub_title = paste("domains: ", domain_cnt);

		print(sub_title);
	
		if(bayes)
		{
			quartz();
			hist_tmp = hist(as.matrix(all$llr), col=color1, main=title, sub=sub_title, xlab="log likelihood ratio", freq=TRUE, breaks=break_points, ylim=ylim, cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2);
			hist_corr_tmp = hist(as.matrix(correct$llr), col=color2, add=TRUE, freq=TRUE, breaks=hist_tmp$breaks);

			ret_val$all= hist_tmp;
			ret_val$correct = hist_corr_tmp;

			#kdrew: remove NaNs from labels because they're stupid looking, replace them with blanks
			labels=hist_corr_tmp$counts/hist_tmp$counts;
			label_nans = is.nan(labels)
			labels=format(labels,digits=2)
			labels[label_nans] = ""
			#labels[labels%in%"NaN"] = ""

			text(hist_tmp$mids, hist_tmp$count,labels=labels,pos=3, cex=2);
		}

		lapply(yield_points, function(yp)
		{
			domain_yield = getYield(eval_table, seq_table, predictors, base_max, base_min, yp);
			abline(v=yp,lty=3);
			if(is.null(ylim))
			{ height = max(hist_tmp$counts); }
			else
			{ height = max(ylim); }
			percent_yield = format(100*(domain_yield/domain_cnt),digits=2);
			text(yp, height, paste(percent_yield,"% domain yield"), pos=4);
		});

		#ret_val$correct = correct;
		#ret_val$all = all;
		ret_val$domain_cnt = domain_cnt;
		return(ret_val);
}

getDomainCount <- function(eval_table, seq_table)
{
		domain_qry = paste("select count(distinct st.seq_id) from ", eval_table, " as et, ",seq_table, " as st where et.domain_sequence_key = st.seq_id");
		print(domain_qry);
		domain_cnt = hpf_qry(domain_qry)[1,1]

		return(domain_cnt);
}


getYield <- function(eval_table, seq_table, predictors, base_max, base_min, yp, yp2=NULL)
{
	domain_qry = paste("select count(distinct st.seq_id) from ", eval_table, " as et, ",seq_table, " as st where et.domain_sequence_key = st.seq_id ");
	domain_qry = paste(domain_qry," and et.",predictors,"_llr >=", yp, " and et.base_llr <= ",base_max," and et.base_llr >= ", base_min, sep="");
	if(!is.null(yp2))
	{
		domain_qry = paste(domain_qry," and et.",predictors,"_llr <", yp2, sep="");
	}
		
	print(domain_qry);
	domain_yield = hpf_qry(domain_qry)[1,1];

	return(domain_yield);
}


funcPredDomain_plotHist <- function(seq_table, eval_table, predictors="s", break_points="Sturges", ylim=NULL, predictor_min = 0, base_max = 0, base_min = -99, color1="grey", color2="seagreen", title="", add=FALSE)
{
	#select sequence_key, MAX(s_llr) as s_llr  from swissprot_mf_seqs_sinceSolved as ss, hpf_results.swissprot_bayes_results as sr where sr.sequence_key = ss.seq_id and s_llr > 0  and base_llr < 0 group by sequence_key
	ret_val = list();
	query = paste(" select MAX(ery.",predictors,"_llr) as llr", sep="");
	query = paste(query, " from ",eval_table, " as ery ");
	query = paste(query, " , ", seq_table," as seq_table");

	query = paste(query, " where seq_table.seq_id = ery.domain_sequence_key");

	query = paste(query, " and ery.base_llr <= ",base_max, " and ery.base_llr >= ",base_min);
	query = paste(query, " and ery.",predictors,"_llr >= ", predictor_min, sep="");
	query = paste(query, " group by ery.domain_sequence_key");

	print(query);
	domain_max_llrs = hpf_qry(query)

	if(title == "")
	{
		title = paste("Histogram of Function Prediction Domain Yield for ", eval_table, ": ",predictors, " predictors"); 
	}

	domain_cnt = getDomainCount(eval_table,seq_table);

	sub_title = paste("domains: ", domain_cnt);

	print(sub_title);
	
	if(add)
	{
		hist_tmp = hist(as.matrix(domain_max_llrs$llr), col=color1, main=title, sub=sub_title, xlab="log likelihood ratio",breaks=break_points, ylim=ylim, add=add);
	}
	else
	{
		quartz();
		hist_tmp = hist(as.matrix(domain_max_llrs$llr), col=color1, main=title, sub=sub_title, xlab="log likelihood ratio",breaks=break_points, ylim=ylim);
	}

	ret_val$domain_max_llrs = hist_tmp;

	#kdrew: remove NaNs from labels because they're stupid looking, replace them with blanks
	labels=hist_tmp$counts;
	label_nans = is.nan(labels)
	labels=format(labels,digits=2)
	labels[label_nans] = ""

	text(hist_tmp$mids, hist_tmp$count,labels=labels,pos=3);

	ret_val$domain_cnt = domain_cnt;
	return(ret_val);
}

funcPredSpecificity_plotHist <- function(seq_table, eval_table, predictors="s", break_points="Sturges", ylim=NULL, predictor_min = 0, base_max = 0, base_min = -99, color1="grey", color2="seagreen", title="", add=FALSE)
{
	#select MIN(ery.base_llr) as llr  from  hpf_results.swissprot_bayes_results_mf  as ery   ,  swissprot_benchmark.swissprot_mf_seqs_sinceSolved  as seq_table  where seq_table.seq_id = ery.sequence_key  and ery.base_llr <=  0  and ery.base_llr >=  -99 and ery.ps_llr >= 2  group by ery.sequence_key

	ret_val = list();
	query = paste(" select MIN(ery.base_llr) as llr", sep="");
	query = paste(query, " from ",eval_table, " as ery ");
	query = paste(query, " , ", seq_table," as seq_table");

	query = paste(query, " where seq_table.seq_id = ery.domain_sequence_key");

	query = paste(query, " and ery.base_llr <= ",base_max, " and ery.base_llr >= ",base_min);
	query = paste(query, " and ery.",predictors,"_llr >= ", predictor_min, sep="");
	query = paste(query, " group by ery.domain_sequence_key");
	print(query);
	base_llrs = hpf_qry(query)

	if(title == "")
	{
		title = paste("Histogram of Function Prediction Specificity for ", eval_table, ": ",predictors, " predictors"); 
	}

	domain_cnt = getDomainCount(eval_table,seq_table);

	sub_title = paste("domains: ", domain_cnt);

	print(sub_title);
	
	if(add)
	{
		hist_tmp = hist(as.matrix(base_llrs$llr), col=color1, main=title, sub=sub_title, xlab="base log likelihood ratio",breaks=break_points, ylim=ylim, add=add);
	}
	else
	{
		quartz();
		hist_tmp = hist(as.matrix(base_llrs$llr), col=color1, main=title, sub=sub_title, xlab="base log likelihood ratio",breaks=break_points, ylim=ylim);
	}

	ret_val$domain_max_llrs = hist_tmp;

	#kdrew: remove NaNs from labels because they're stupid looking, replace them with blanks
	labels=hist_tmp$counts;
	label_nans = is.nan(labels)
	labels=format(labels,digits=2)
	labels[label_nans] = ""

	text(hist_tmp$mids, hist_tmp$count,labels=labels,pos=3);

	ret_val$domain_cnt = domain_cnt;
	return(ret_val);
}

funcPredCompare <-function(domainYield1, domainYield2, title="", subtitle="", color1="seagreen", color2="orange")
{
	quartz();
	plot(domainYield1$domain_max_llrs$mids, domainYield1$domain_max_llrs$counts, type="b", col=color1, main=title, sub=subtitle, xlab="log likelihood ratio", ylab="number of domains", cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2, pch=19);
	points(domainYield2$domain_max_llrs$mids, domainYield2$domain_max_llrs$counts, type="b", col=color2, cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2, pch=23);
	#text(domainYield1$domain_max_llrs$mids, domainYield1$domain_max_llrs$counts, format(domainYield1$domain_max_llrs$counts/domainYield1$domain_cnt,digits=2), pos=3);
	#text(domainYield2$domain_max_llrs$mids, domainYield2$domain_max_llrs$counts, format(domainYield2$domain_max_llrs$counts/domainYield2$domain_cnt,digits=2), pos=3);
	text(domainYield1$domain_max_llrs$mids, domainYield1$domain_max_llrs$counts, domainYield1$domain_max_llrs$counts, pos=3, cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2);
	text(domainYield2$domain_max_llrs$mids, domainYield2$domain_max_llrs$counts, domainYield2$domain_max_llrs$counts, pos=3, cex=2, cex.lab=2, cex.main=2, cex.sub=1.5, cex.axis=2);
}
