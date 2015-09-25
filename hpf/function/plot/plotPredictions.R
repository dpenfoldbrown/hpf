# Creating a string concatenation operator
"||" <- function(...) UseMethod("||")
"||.default" <- .Primitive("||")
"||.character" <- function(...) paste(...,sep="")

source("sql_setup.R");
source("go_hpf_utils.R");
#source("evaluate_predictions.R");

# expression(paste(x,symbol("\307"),y)
symbol_intersect = "\307"
symbol_union = "\310"
symbol_complement = "\\"

plot_function = function(experiment="euk",random=3000,unannotated=T){
	#op <- par(mfrow = c(2, 2),pty = "s")
	colours = list("red","green","black","brown")
	legend_names = c("PC","PCS","random")

#kdrew: things in comments are just so the script does not take a long time to run, in the time I wrote this comment I could have created an if block
	
#	domain_type=c("psiblast")
#        #kdrew: removed the medium_spec() call
#	data = load(domain_type_seqs(protein_seqs(experiment=experiment,with_mf=T,with_scop=T,random=random,domain_type=domain_type),with_scop=T,domain_type=domain_type), unannotated=T,base_max=-4)
#	randdata = load(domain_type_seqs(protein_seqs(experiment=experiment,with_mf=T,with_scop=T,random=random,domain_type=domain_type),with_scop=T,domain_type=domain_type), unannotated=T,base_max=-4,random=T)
#	perf = ir_graph(data,randdata,plot=F,names=c("pl_llr","pls_llr"))
#	limit = perf@x.values[[2]][max(which(perf@alpha.values[[2]] > -99))]
#	pdf(experiment||"_prec_rec_psi.pdf")
#	plot(perf, col=colours,lwd=2, xlim=c(0,limit))
#	title("Eukaryote PDB-Blast Domains")
#	#legend(0.5,1.0,legend_names,col=unlist(colours),lwd=2)
#	dev.off()
	
	domain_type=c("fold_recognition")
	#kdrew: removed the medium_spec() call
	data = load(domain_type_seqs(protein_seqs(experiment=experiment,with_mf=T,with_scop=T,random=random,domain_type=domain_type),with_scop=T,domain_type=domain_type), unannotated=T,base_max=-4)
	randdata = load(domain_type_seqs(protein_seqs(experiment=experiment,with_mf=T,with_scop=T,random=random,domain_type=domain_type),with_scop=T,domain_type=domain_type), unannotated=T,base_max=-4,random=T)
	perf = ir_graph(data,randdata,plot=F,names=c("pl_llr","pls_llr"))
	limit = perf@x.values[[2]][max(which(perf@alpha.values[[2]] > -99))]
	pdf(experiment||"_prec_rec_fr.pdf")
	plot(perf, col=colours,lwd=2,xlim=c(0,limit))
	title("Eukaryote Fold Recognition Domains")
	#legend(0.5,1.0,legend_names,col=unlist(colours),lwd=2)
	dev.off()	

#	domain_type=c("msa","pfam","unassigned")
	#kdrew: removed the medium_spec() call
#	data = load(domain_type_seqs(protein_seqs(experiment=experiment,with_mf=T,with_scop=T,random=random,domain_type=domain_type),with_scop=T,domain_type=domain_type), unannotated=T,base_max=-4)
#	perf = ir_graph(data,plot=F,names=c("pl_llr","pls_llr","pls_random"))
#	limit = perf@x.values[[2]][max(which(perf@alpha.values[[2]] > -99))]
#	pdf(experiment||"_prec_rec_dn.pdf")
#	plot(perf, col=colours,lwd=2,xlim=c(0,limit))
#	title("Eukaryote De Novo Domains")
#	legend(0.5,0.9,legend_names,col=unlist(colours),lwd=2)
#	dev.off()
#	#par(op)

#	domain_type=c("msa","pfam","unassigned","fold_recognition","psiblast")
#	#domain_type=c("unassigned")
#	#kdrew: removed the medium_spec() call
#	data = load(domain_type_seqs(protein_seqs(experiment="sap",with_mf=T,with_scop=T,random=random,domain_type=domain_type,go_table="since_solved.since_solved_iea_golite_062009"),with_scop=T,domain_type=domain_type), unannotated=T, bayes_table="since_solved.bayes_golite_062009_3_redux",go_table="since_solved.since_solved_iea_golite_062009",base_max=-4 )
#	randdata = load(domain_type_seqs(protein_seqs(experiment="sap",with_mf=T,with_scop=T,random=random,domain_type=domain_type,go_table="since_solved.since_solved_iea_golite_062009"),with_scop=T,domain_type=domain_type), unannotated=T, bayes_table="since_solved.bayes_golite_062009_3_redux",go_table="since_solved.since_solved_iea_golite_062009",base_max=-4, random=T )
#	perf = ir_graph(data,randdata, plot=F,names=c("pl_llr","pls_llr"))
#	limit = perf@x.values[[2]][max(which(perf@alpha.values[[2]] > -99))]
#	pdf("sap_prec_rec_dn.pdf")
#	plot(perf, col=colours,lwd=2,xlim=c(0,limit))
#	title("SAP De Novo Domains")
#	legend(0.5,0.9,legend_names,col=unlist(colours),lwd=2)
#	dev.off()
}

plot_orthogonal = function(data, type1="rec", type2="prec",rate=0.5,rand=500){
	#op <- par(mfrow = c(2, 2),pty = "s")
	xlim=c(1.0,0.5)
	o = orthogonal_predictions(data,predictor1="pl_llr",predictor2="pls_llr",set="mf_acc",type1=type1,type2=type2,rate=rate,plot=F,rand=rand)
	ylim = c(0,max(o[[2]]+o[[3]]+o[[4]]))
	pdf("orth_mf.pdf")
	plot(o[[1]],o[[2]]+o[[3]],xlab="Precision",ylab="Unique Molecular Functions",col="black",type='l',xlim=xlim, lwd=2, ylim=ylim)
	lines(o[[1]],o[[3]],col="red", lwd=2)
	lines(o[[1]],o[[2]]+o[[3]]+o[[4]],col="blue", lwd=2)
	#legend(1,ylim[2],
	#		expression(paste("PC ",symbol("\307")," PCS"),"PC \\ PCS","PCS \\ PC"),
	#		col=c("black","red","blue"),lwd=2)
	dev.off()
	
	#o = orthogonal_predictions(data,predictor1="pl_llr",predictor2="s_llr",set="mf_acc",type1=type1,type2=type2,rate=rate,plot=F)
	#ylim = c(0,max(o[[2]],o[[3]],o[[4]]))
	#plot(o[[1]],o[[2]],xlab="Precision",ylab="Unique Molecular Functions",col="black",type='l',xlim=xlim, lwd=2, ylim=ylim)
	#lines(o[[1]],o[[3]],col="red", lwd=2)
	#lines(o[[1]],o[[4]],col="blue", lwd=2)
	#legend(1,ylim[2],
	#		expression(paste("PL ",symbol("\307")," S"),"PL \\ S","S \\ PL"),
	#		col=c("black","red","blue"),lwd=2)
	
	o = orthogonal_predictions(data,predictor1="pl_llr",predictor2="pls_llr",set="parent_sequence_key",type1=type1,type2=type2,rate=rate,plot=F,rand=rand)
	ylim = c(0,max(o[[2]]+o[[3]]+o[[4]]))
	pdf("orth_seq.pdf")
	plot(o[[1]],o[[2]]+o[[3]],xlab="Precision",ylab="Unique Proteins",col="black",type='l',xlim=xlim, lwd=2, ylim=ylim)
	lines(o[[1]],o[[3]],col="red", lwd=2)
	lines(o[[1]],o[[2]]+o[[3]]+o[[4]],col="blue", lwd=2)
	#legend(1,ylim[2],
	#		expression(paste("PC ",symbol("\307")," PCS"),"PC \\ PCS","PCS \\ PC"),
	#		col=c("black","red","blue"),lwd=2)
	dev.off()
	
	#o = orthogonal_predictions(data,predictor1="pl_llr",predictor2="s_llr",set="parent_sequence_key",type1=type1,type2=type2,rate=rate,plot=F)
	#ylim = c(0,max(o[[2]],o[[3]],o[[4]]))
	#plot(o[[1]],o[[2]],xlab="Precision",ylab="Unique Proteins",col="black",type='l',xlim=xlim, lwd=2,ylim=ylim)
	#lines(o[[1]],o[[3]],col="red", lwd=2)
	#lines(o[[1]],o[[4]],col="blue", lwd=2)
	#legend(1,ylim[2],
	#		expression(paste("PL ",symbol("\307")," S"),"PL \\ S","PLS \\ PL"),
	#		col=c("black","red","blue"),lwd=2)
	
	#par(op)
}


label = function(data){
	r = split(data)
	return(rbind(r[[1]],r[[2]],r[[3]]))
}


#' 
#' @param data 
#' @param predictor1 
#' @param predictor2 
#' @returnType 
#' @return 
#' @author patrick
#' @export
orthogonal_predictions = function(data, predictor1="pl_llr", predictor2="pls_llr",rand=500, plot=T,y=NA,x=0.1,set="mf_acc",rate=0.5,type1="prec",type2="rec"){
	a = ir_graph(data,names=c(predictor1,predictor2),plot=F,type1=type1,type2=type2)
	recall_primary = a@x.values[[1]]
	recall_other = a@x.values[[2]]
	alpha_primary = a@alpha.values[[1]]
	alpha_other = a@alpha.values[[2]]
	diff1 = c()
	diff2 = c()
	union_count = c()
	intersection_count = c()
	threshold = recall_primary[which(recall_primary >= rate)]
	s = sort(sample(threshold,min(rand,length(threshold))), decreasing=T)
	#print("threshold"||threshold)
	#print("rand"||rand)
	i = length(s)
	print("length"||i)
	for(recall in s){
		thresh_primary = min(alpha_primary[which(recall_primary >= recall)])
		mf_primary = unique(data[which(data[[predictor1]] >= thresh_primary), ][[set]])
		thresh_other = min(alpha_other[which(recall_other >= recall)])
		mf_other = unique(data[which(data[[predictor2]] >= thresh_other), ][[set]])
		intersection = intersect(mf_primary, mf_other)
		uni = union(mf_primary, mf_other)
		diff = setdiff(uni, intersection)
		
		diff1 = c(diff1, length(setdiff(mf_primary,mf_other)))
		diff2 = c(diff2, length(setdiff(mf_other,mf_primary)))
		union_count = c(union_count,length(uni))
		intersection_count = c(intersection_count,length(intersection))
		i = i-1
		if(i %% 50 == 0){
			print(""||i||" "||predictor1||":"||length(setdiff(mf_primary,mf_other))||" "||predictor2||":"||length(setdiff(mf_other,mf_primary))||" intersection:"||length(intersection)||" precision:"||recall)
		}
	}
	if(plot){
		if(is.na(y)){
			y = max(intersection_count)
		}
		plot(s,intersection_count,xlab="Recall",ylab="Size of "||set||" set",col="black",type='l',xlim=c(max(s),min(s)))
		lines(s,diff1,col="red")
		lines(s,diff2,col="blue")
		legend(x,y,c("intersection("||predictor1||","||predictor2||")",predictor1||" - "||predictor2,predictor2||" - "||predictor1),col=c("black","red","blue"),lwd=2)
	}
	return(list(s,intersection_count,diff1,diff2))
}

genome_model = function(domain_types = list("psiblast","fold_recognition",c("msa","pfam","unassigned")),domain_table=NA,experiment=NA,random=NA,predictor="pls_llr",min_accuracy=0.5,spec_bins=c(1,2,3),bayes_table="functionPredictions.bayes_golite_062009_3",go_table="hddb_IEA.hddb_IEA_goLite_062009",distance=0){
	type_thresholds = list()
	for(domains in domain_types){
		print("Generating Model for "||paste(domains,collapse=","))
		domain_mins = type_mins(domain_type=domains,domain_table=domain_table, experiment=experiment,random=random,predictor=predictor,min_accuracy=min_accuracy,spec_bins=spec_bins,bayes_table=bayes_table,go_table=go_table,distance=distance)
		type_thresholds = c(type_thresholds,list(domain_mins))
	}
	return(type_thresholds)
}

genome_coverage = function(error_model,data=NA,experiment=NA,random=NA,with_mf=NA,with_scop=F,predictor="pls_llr",spec_bins=c(2,3),domain_types = list("psiblast","fold_recognition",c("msa","pfam","unassigned")),bayes_table="functionPredictions.bayes_golite_062009_3",domain_table=NA,go_table="hddb_IEA.hddb_IEA_goLite_062009"){
	
	if(is.na(data)){
		print("Loading global samples for "||paste(experiment,collapse=","))
		if(is.na(domain_table)){
			domain_table = domain_type_seqs(protein_seqs(experiment=experiment,random=random,with_mf=with_mf,with_scop=with_scop,go_table=go_table))
		}
		data = label(load(domain_table,bayes_table=bayes_table,go_table=go_table))
	}
	print("Found "||length(unique(data$parent_sequence_key))||" proteins in data")
	bindata = bin_llr(data)
	keep = data.frame()
	for(i in 1:length(domain_types)){
		for(bin in spec_bins){
			domains = domain_types[[i]]
			domain_min = error_model[[i]][bin]
			print("Selecting "||paste(domains,collapse=","))
			domain_data = bindata[[bin]][which(bindata[[bin]]$domain_type %in% domains),]
			print("  "||length(unique(domain_data$parent_sequence_key))||" proteins, "||length(unique(domain_data$domain_sequence_key))||" domains, "||length(domain_data[[predictor]])||" predictions")
			print("Thresholding "||predictor||">="||domain_min)
			domain_keep = domain_data[which(domain_data[[predictor]] >= domain_min),]
			print("  Keeping "||length(unique(domain_keep$parent_sequence_key))||" proteins")
			keep = rbind(keep,domain_keep)
		}
	}
	venn_type(keep)
	return(list(data,keep))
}

type_mins = function(domain_type=NA,experiment=NA,random=NA,predictor="pls_llr",min_accuracy=0.5,spec_bins=c(1,2,3),data=NA,domain_table=NA,bayes_table="functionPredictions.bayes_golite_062009_3",go_table="hddb_IEA.hddb_IEA_goLite_062009",distance=0){
	#domain_types = list("fold_recognition","psiblast",c("msa","pfam","unassigned"))
	if(!is.na(domain_table)){
		print("DOMAIN TABLE")
		with_mf = label(load(domain_table,bayes_table=bayes_table,go_table=go_table,distance=distance))
	} else if(is.na(data)){
		print("Loading Data")
		domain_table = domain_type_seqs(protein_seqs(experiment=experiment,domain_type=domain_type,with_mf=T,random=random,with_scop=T,go_table=go_table),domain_type=domain_type)
		print("Domain Table: "||domain_table)
		with_mf = load(domain_table,bayes_table=bayes_table,unannotated=F,go_table=go_table,distance=distance)
		#print(length(with_mf$pls_llr))
		with_mf = label(with_mf)
		#print(length(with_mf$pls_llr))
	} else {
	       print("existing data")
		with_mf = label(data)
	}
	# Split data up by specificity, model medium to specific functions
	bins = bin_llr(with_mf)
	mins = c()
	for(bin in spec_bins){
		d = bins[[bin]]
		print(unique(d$labels))
		print(length(d$labels))
		predictor_min = precision_threshold(d,predictor,min_accuracy=min_accuracy)
		mins = c(mins,predictor_min)
	}
	return(mins)
}

#' Generates precision threshold graphs for all domain types and experiments.
#' @param med 
#' @returnType 
#' @return 
#' @author Patrick
#' @export
genome_prec_graph = function(med=T){
	i = 0
	colours = list("purple","red","orange","blue","black","brown")
	legend_names = c()
	legend_cols = c()
	for(exp in list("prok","euk")){
		for(domain_type in list("psiblast","fold_recognition",c("msa","pfam","unassigned"))){
			i = i+1
			col = colours[[i]]
			if(i > 1){
				par(new=T)
				axes=F
				xlab=NA
				ylab=NA
			} else {
				axes=T
				xlab="pls_llr"
				ylab="precision"
			}
			data = load(domain_type_seqs(protein_seqs(experiment=exp,with_scop=T,with_mf=T,random=5000,domain_type=domain_type),domain_type=domain_type))
			if(med){
				data = medium_spec(data)
			}
			a = ir_graph(data,"prec","rec",plot=F)
			plot(a@alpha.values[[7]],a@y.values[[7]],
					xlab=xlab,ylab=ylab,
					type="l",col=col,xlim=c(-5,10),ylim=c(0,1))
			if(domain_type==c("msa","pfam","unassigned")){
				domain_type="mcm"
			}
			legend_names=c(legend_names,exp||"_"||domain_type)
			legend_cols=c(legend_cols,col)
		}		
	}
	legend(-5,1,legend_names,col=legend_cols,lwd=2)
	if(med){
		spec=" Specific"
	} else {
		spec=""
	}
	title("Precision/Threshold"||spec||" Functions")
}

sccs_graph = function(data,sf,location="/home/patrick/Sites/hpf/plot/mcm/sf/"){#,bayes_table="functionPredictions.bayes_golite_062009_3"){
		print(sf)
		#domain = "(select distinct d.domain_sequence_key,d.domain_type from hpf.domain_sccs d join hddb_IEA.hddb_iea_golite_062009 g on d.parent_sequence_key=g.sequence_key and g.term_type='molecular_function' and g.acc!='GO:0003674' where d.domain_type in ('msa','pfam','unassigned') and substring_index(d.sccs,'.',3)='"||sf||"')"
		#data = label(load(domain,bayes_table=bayes_table))
		png(location||sf||".png",bg="white",width=800,height=800)
		ir_graph(data)
		title(""||sf||" All Functions")
		dev.off()
		png(location||sf||".thresh.png",bg="white",width=800,height=800)
		prec_graph(data)
		title(""||sf||" Precision vs 'pls_llr' All Functions")
		dev.off()
		png(location||sf||".spec.png",bg="white",width=800,height=800)
		ir_graph(medium_spec(data))
		title(""||sf||" Specific Functions")
		dev.off()
		
	
}

# mcm_euk_data2 = load(domain_type_seqs(protein_seqs(experiment="euk",with_scop=T, domain_type=c("msa","unassigned","pfam")),domain_type=c("msa","unassigned","pfam"),with_scop=T))
# med_spec_data2 = medium_spec(mcm_euk_data2)
# dsk_sampled100 = sample(med_spec_data2$domain_sequence_key, 100)
# med_spec_data2_100 = med_spec_data2[med_spec_data2$domain_sequence_key %in% dsk_sampled100,]
# bar_graph(med_spec_data2_100,paper_graph=T)

# fr_euk_data2 = load(domain_type_seqs(protein_seqs(experiment="euk",with_scop=T, domain_type=c("fold_recognition")),domain_type=c("fold_recognition"),with_scop=T))
# med_spec_fr_data2 = medium_spec(fr_euk_data2)

# fr_dsk_sampled100 = sample(med_spec_fr_data2$domain_sequence_key, 100)
# med_spec_fr_data2_100 = med_spec_fr_data2[med_spec_fr_data2$domain_sequence_key %in% fr_dsk_sampled100,]
# bar_graph(med_spec_fr_data2_100,paper_graph=T)

bar_graph = function(data, sampled=100, ordered_by = "pl_llr", threshold = -99){
	require(ROCR)

	#dsk_sampled = sample(data$domain_sequence_key, sampled)
	#data_sampled = data[data$domain_sequence_key %in% dsk_sampled,]

	thres_data = data[data$pls_llr > threshold,]

	row_sampled = sample(seq(nrow(thres_data)),sampled)
	data_sampled = thres_data[row_sampled,]
	

	data_l = label(data_sampled)
	data_l$pls_random = sample(data_l$pls_llr)

	compare_pl_llr = data_l[do.call(order, -data_l[ordered_by]),]$pl_llr
	compare_pls_llr = data_l[do.call(order, -data_l[ordered_by]),]$pls_llr

	par(mfrow=c(2,1))
	plot_bar_graph(compare_pl_llr)
	title(main="Molecular Function predictions")
	title(ylab="pl_llr")
	plot_bar_graph(compare_pls_llr)
	title(ylab="pls_llr")
	title(xlab=paste("index ordered by: ",ordered_by))
}

plot_bar_graph = function(compare_llr)
{
	w = .5
	y = 1:length(compare_llr)

	plot.new()
	plot.window(xlim=c(-1, length(compare_llr)+5), ylim=c(-9, 10))
	yticks <- seq(-9,10, 3)
	xticks <- seq(1,length(compare_llr)+5, 200)
	axis(2, at=yticks, labels=yticks, pos=0)
	axis(1, at=xticks, labels=xticks, pos=-9)

	rect(y-w,0,y+w,compare_llr)

}

ir_graph = function(data, randdata=NULL, type1="prec",type2="rec",x=0.7,y=1,plot=T, names = c("p_llr","l_llr","s_llr","ps_llr","ls_llr","pl_llr","pls_llr","pls_random"), colours = list("yellow","purple","red","blue","green","orange","black","brown"),paper_graph=F){
	require(ROCR)
	if(paper_graph){
		names = c("s_llr","pl_llr","pls_llr","pls_random")
		legend_names = c("structure", "process/component", "process/component/structure", "random")
		colours = list("red","green","black","brown")
	} else {
		legend_names = names
	}
	data_l = label(data)
	print(data_l)
	#kdrew: this is not correct because it only samples scores above the pls_llr cutoff and not the full prediction score
	data_l$pls_random = sample(data_l$pls_llr)
	#data_l$pls_random = data_l$base_llr
	print(data_l)
	pp = list()
	ll = list()
	for(col in names){
	#	print("col"||col)
	#	print(list(data_l[[col]]))
	#	print(list(data_l$labels))
		pp = c(pp, list(data_l[[col]]))
		ll = c(ll, list(data_l$labels))
	}
	if(!is.null(randdata))
	{
		data_r = label(randdata)
		print("random data")
		print(data_r)
		print(length(data_r[["base_llr"]]))
		print(length(data_r$labels))
		pp = c(pp,list(data_r[["base_llr"]]))
		ll = c(ll,list(data_r$labels))
	}
	pred = prediction(pp,ll)
	#print(pred)
	perf = performance(pred,type1,type2)
	#perf = performance(pred,"tpr","fpr")
	if(plot){
		plot(perf, col=colours)
		legend(x,y,legend_names,col=unlist(colours),lwd=1)
	}
	#kdrew: kinda meaningless because it does not work for prec/recall curves only ROC
	auc_perf = performance(pred,'auc')
	print("AUC:"|| auc_perf@y.values)
	return(perf)
}

#' Processes all of the superfamilies that have known molecular function.
#' Calculates prec/rec values for PL and PL with Structure to identify
#' superfamilies with structural improvement.
#' @param sccs_query 
#' @returnType 
#' @return 
#' @author patrick
#' @export
prec_sccs = function(min_accuracy=0.5,sccs_query="select distinct sccs,num,avg_conf from (select count(distinct d.domain_sequence_key) as num, substring_index(sccs,'.',3) as sccs, avg(confidence) as avg_conf from hpf.protein p join domain_sccs d on p.sequence_key=d.parent_sequence_key join hddb_IEA.hddb_iea_golite_062009 g on d.parent_sequence_key=g.sequence_key and g.term_type='molecular_function' and g.acc != 'GO:0003674' where domain_type in ('msa','pfam','unassigned') and p.experiment_key not in (849,824,890,901,6,18,887) group by substring_index(sccs,'.',3)) p where p.num > 10",bayes_table="functionPredictions.bayes_golite_062009_3"){
	print(sccs_query)
	sccs_data = hpf_qry(sccs_query)
	print("Number of superfamilies to process "||length(sccs_data[["sccs"]]))
	#return(sccs_data)
	for(i in 1:length(sccs_data[["sccs"]])){
		print("i "||i)
		sccs = sccs_data[i, ][["sccs"]]
		num = sccs_data[i, ][["num"]]
		avg_conf = sccs_data[i, ][["avg_conf"]]
		#print(sccs)
		#print(num)
		#print(avg_conf)
		print(sccs||" "||num||" "||avg_conf)
		
		domains = "(select distinct d.domain_sequence_key,d.domain_type from hpf.protein p join hpf.domain_sccs d on p.sequence_key=d.parent_sequence_key join hddb_IEA.hddb_iea_golite_062009 g on d.parent_sequence_key=g.sequence_key and g.term_type='molecular_function' and g.acc!='GO:0003674' where d.domain_type in ('msa','pfam','unassigned') and p.experiment_key not in (849,824,890,901,6,18,887) and substring_index(d.sccs,'.',3)='"||sccs||"')"
		data = label(load(domains,bayes_table=bayes_table))
		sdata = medium_spec(data)
		#tryCatch({
		#					sccs_graph(data,sccs)
		#		}, error=function(e) {
		#			print("caught an error making graph")
		#		})
		
		p = precision_rec(data,predictor="p_llr",min_accuracy=min_accuracy)
		s = precision_rec(data,predictor="s_llr",min_accuracy=min_accuracy)
		ps = precision_rec(data,predictor="ps_llr",min_accuracy=min_accuracy)
		
		pl = precision_rec(data,predictor="pl_llr",min_accuracy=min_accuracy)
		pls = precision_rec(data,predictor="pls_llr",min_accuracy=min_accuracy)
		
		s_pl = precision_rec(sdata,predictor="pl_llr",min_accuracy=min_accuracy)
		s_pls = precision_rec(sdata,predictor="pls_llr",min_accuracy=min_accuracy)
		query = "insert ignore into sf_func_rec (sccs,domains,avg_conf,prec,pl,pls,s_pl,s_pls,bayes_table,p,s,ps) values ('"||sccs||"',"||num||","||avg_conf||","||min_accuracy||","||pl||","||pls||","||s_pl||","||s_pls||",'"||bayes_table||"',"||p||","||s||","||ps||")"
		print(query)
		hpf_qry(query)
	}
}

prec_graph = function(data, predictor="pls_llr",col="black"){
	data = label(data)
	pred = prediction(data[[predictor]],data$labels)
	perf = performance(pred,"prec","rec")
	plot(perf@alpha.values[[1]],perf@y.values[[1]],col=col,type="l",xlim=c(-3,10),ylim=c(0,1.0))
}
	

medium_spec = function(data){
	bins = bin_llr(data)
	return(rbind(bins[[2]],bins[[3]]))
}


protein_stat = function(data,terms){
	data = label(data) # remove non-predictions
	stat = data.frame()
	for(protein_key in unique(data$parent_sequence_key)){
		#print(protein_key)
		p_data = data[which(data$parent_sequence_key %in% protein_key), ]
		t = terms[which(terms$sequence_key %in% protein_key), ]
		
		#domains = length(unique(p_data$domain_sequence_key))
		query = "select sum(size) as domains from domain_sccs where parent_sequence_key="||protein_key
		domains = hpf_qry(query)$domains
		
		predictions = length(p_data$mf_acc)
		unique_mf = length(unique(p_data$mf_acc))
		t_num = t$num_terms
		if(length(t$num_terms)==0){
			t_num=0
		}
		
		#print("domains"||domains)
		#print("predictions"||predictions)
		#print("unique"||unique_mf)
		#print("terms"||t_num)
		
		p_stat = data.frame(parent_sequence_key=protein_key,
					domains=domains,
					predictions=predictions,
					unique_mf=unique_mf,
					terms=t_num)
		#print(p_stat)
		stat = rbind(stat,p_stat)
	}
	return(stat)
}

protein_terms = function(protein_table="hpf_test.protein_seqs", terms=c("biological_process","cellular_component"),sequence_key="parent_sequence_key",go_table="hddb_IEA.hddb_iea_golite_062009"){
	query = 
		"select "||sequence_key||" as sequence_key, count(distinct t1.acc) as num_terms
			from "||protein_table||" seq_table
			left outer join "||go_table||" g
				on seq_table."||sequence_key||"=g.sequence_key 
				and g.evidence_code in ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA','hpf_IEA')
				and g.acc != 'GO:0003674'
				and g.acc != 'GO:0005575'
				and g.acc != 'GO:0008150'
			left outer join mygoLite_062009.term t2
				on g.acc=t2.acc 
				and t2.term_type in ('"||paste(terms,collapse="','")||"')
				and t2.is_obsolete = 0
			left outer join mygoLite_062009.graph_path path
				on t2.id = path.term2_id
			left outer join mygoLite_062009.term t1
				on path.term1_id=t1.id
				and t1.is_obsolete = 0
			group by seq_table."||sequence_key

	
	
	print(query)
	res = hpf_qry(query)
	return(res)
}

protein_terms_metric = function(protein_table="hpf_test.protein_seqs", terms=c("biological_process","cellular_component"),sequence_key="parent_sequence_key",go_table="hddb_IEA.hddb_iea_golite_062009", functionTable="functionTables.log_ratio_goLite_062009_3"){
	query = 
		"select distinct a.sequence_key,a.acc,ft.metric from 
			(select distinct "||sequence_key||" as sequence_key, t1.acc
			from "||protein_table||" seq_table
			join "||go_table||" g
				on seq_table."||sequence_key||"=g.sequence_key 
				and g.evidence_code in ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA','hpf_IEA')
				and g.acc != 'GO:0003674'
				and g.acc != 'GO:0005575'
				and g.acc != 'GO:0008150'
			join mygoLite_062009.term t2
				on g.acc=t2.acc 
				and t2.term_type in ('"||paste(terms,collapse="','")||"')
				and t2.is_obsolete = 0
			join mygoLite_062009.graph_path path
				on t2.id = path.term2_id
			join mygoLite_062009.term t1
				on path.term1_id=t1.id
				and t1.is_obsolete = 0 ) a
			join "||functionTable||" ft
			     	on a.acc=ft.acc and ft.acc2 is NULL"
	
	print(query)
	res = hpf_qry(query)
	return(res)
}


domain_types = function(data){
	data_psi = data[which(data$domain_type %in% "psiblast"), ]
	data_fr = data[which(data$domain_type %in% "fold_recognition"), ]
	data_mcm = data[which(data$domain_type %in% c("msa","pfam","unassigned")), ]
	return(list(data_psi,data_fr,data_mcm))
}

precision_threshold = function(data, predictor="pls_llr", min_accuracy=0.5){
	a = label(data)
	require(ROCR)
	pred_min = 100
	tryCatch({
		perf = performance(prediction(a[[predictor]],a$labels),"prec","rpp")
		# Bypass the first couple because those are unreliable
		length = length(perf@alpha.values[[1]])
		pred_min = min(perf@alpha.values[[1]][100:length][which(perf@y.values[[1]][100:length] >= min_accuracy)])
	}, error=function(e) {
		print("caught an error")
	})
	return(pred_min)
}

precision_rec = function(data, predictor="pls_llr", min_accuracy=0.5){
	a = label(data)
	require(ROCR)
	rec_min = 0
	tryCatch({
				perf = performance(prediction(a[[predictor]],a$labels),"prec","rec")
				length = length(perf@alpha.values[[1]])
				rec_min = max(perf@x.values[[1]][100:length][which(perf@y.values[[1]][100:length] >= min_accuracy)])
			}, error=function(e) {
				print("caught an error")
			})
	if(rec_min == Inf | rec_min==-Inf){
		rec_min=0
	}
	return(rec_min)
}


venn_type = function(data,type="parent_sequence_key"){
	require(limma)
	universe = unique(data[[type]])
	r = domain_types(data)
	sets = list()
	names = c("psiblast","fold_recognition","mcm")
	for(i in 1:length(r)){
		sets = c(sets,list(unique(r[[i]][[type]])))
	}
	Counts <- matrix(0, nrow=length(universe), ncol=length(names))
	colnames(Counts) <- names
	print("counts")
	for (i in 1:length(universe)){
		for (j in 1:length(sets)){
			Counts[i,j] <- universe[i] %in% sets[[j]]
		}
	}
	vennDiagram( vennCounts(Counts) )
}

venn_predictors = function(data, type="parent_sequence_key", min_accuracy=0.5, col_mins=NA){
	require(limma)
	a = label(data)
	names = c("p_llr","l_llr","pls_llr")
	sets = list()
	universe = c()
	for(i in 1:length(names)){
		col=names[i]
		print("col "||col)
		if(is.na(col_mins)){
			pred_min = precision_threshold(a,col,min_accuracy=min_accuracy)
		} else {
			pred_min = col_mins[i]
		}
		print("min "||pred_min)
		if(pred_min == Inf)
			pred_type=c()
		else
			pred_type = unique(a[which(a[[col]] >= pred_min),][[type]])
		sets = c(sets,list(pred_type))
		universe = sort(union(universe,pred_type))
	}
	Counts <- matrix(0, nrow=length(universe), ncol=length(names))
	colnames(Counts) <- names
	print("counts")
	for (i in 1:length(universe)){
		for (j in 1:length(sets)){
			Counts[i,j] <- universe[i] %in% sets[[j]]
		}
	}
	vennDiagram( vennCounts(Counts) )
}

protein_seqs = function(experiment=NA, with_mf=T, with_scop=F, random=NA, domain_type=NA,go_table="hddb_IEA.hddb_IEA_goLite_062009"){
	if(experiment=="prok"){
		experiment=c(29, 30, 31, 34, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885)#, 889, 920)
	}
	else if(experiment=="euk"){
		experiment=c(804, 822, 823, 825, 826, 827, 886, 915, 924, 900)#888, 917,4
	} else if (experiment=="all"){
		experiment=c(29, 30, 31, 34, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 804, 822, 823, 825, 826, 827, 886, 915, 924, 900)
	}
	
	# Delete the old table
	query = "drop temporary table if exists proteins"
	print(query)
	hpf_qry(query)
	# Build a table of proteins
	domain_table="domain"
	if(with_scop)
		#domain_table="since_solved.p_domain_sccs_redux"
		domain_table="domain_sccs"
	where = F
	query = "select distinct p.sequence_key as parent_sequence_key 
		from "||domain_table||" d join protein p on d.parent_sequence_key=p.sequence_key"
	
	if(!is.na(domain_type)){
		if(!where)
			query = query||" where"
		else
			query = query||" and"
		query = query||" d.domain_type in ('"||paste(domain_type,collapse="','")||"')"
		where = T
	}
	if(with_scop){
		if(!where)
			query = query||" where"
		else
			query = query||" and"
		query = query||" d.sccs is not NULL"
		where = T
	}
	if(!is.na(experiment)){
		if(!where)
			query = query||" where"
		else
			query = query||" and"
		query = query||" p.experiment_key in ("||paste(experiment,collapse=",")||")";
	}
	if(experiment == "sap")
	{
		query = "select distinct sequence_key as parent_sequence_key from since_solved.p_mcm_redux"
	}
	query = "create temporary table proteins "||query
	print(query)
	hpf_qry(query)
	
	query = "create index parent_sequence_key on proteins(parent_sequence_key)"
	print(query) 
	hpf_qry(query)
	
	
	# Now remove tainted protein sequences
	print("Before removing experiments "||hpf_qry("select count(distinct parent_sequence_key) from proteins")||" proteins")
	# Because of the many-many mapping with experiment/protein we need to be
	# be careful to delete these
	query = "delete from proteins where parent_sequence_key in 
			(select p.sequence_key from hpf.protein p where p.experiment_key in (849,824,890,901,6,18,887))"
	print(query)
	hpf_qry(query)
	
	
	# Delete the old table
	query = "drop temporary table if exists protein_mf"
	print(query)
	hpf_qry(query)
	# Create a mapping to their molecular function
	query = "create temporary table protein_mf
		select p.parent_sequence_key, g.* from
		proteins p
		join "||go_table||" g
			on p.parent_sequence_key=g.sequence_key 
			and g.evidence_code in ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA','hpf_IEA')
			and g.acc != 'GO:0003674'
		join mygoLite_062009.term t 
			on g.acc=t.acc 
			and t.term_type='molecular_function'
			and t.is_obsolete = 0"
	print(query)
	hpf_qry(query)
	query = "create index parent_sequence_key on protein_mf(parent_sequence_key)"
	print(query)
	hpf_qry(query)
	
	# Delete the old table
	query = "drop table if exists hpf_test.protein_seqs"
	print(query)
	hpf_qry(query)

	# Now create the final table
	query = "create table hpf_test.protein_seqs select distinct p.parent_sequence_key 
		from proteins p 
		left outer join protein_mf m
			on p.parent_sequence_key=m.parent_sequence_key"
	if(is.na(with_mf)){
		#do nothing
	}
	else if(with_mf==T){
		#print("With mf")
		query = query||" where m.acc is not NULL"
	} else if (with_mf==F) {
		where = T
		query = query ||" where m.acc is NULL"
	}	
	if(!is.na(random)){
		query = query||" order by rand() limit "||random
	}
	print(query)
	hpf_qry(query)
	
	query = "create index parent_sequence_key on hpf_test.protein_seqs (parent_sequence_key)"
	hpf_qry(query)
	#query = "create index (acc) on hpf_test.protein_seqs (acc)"
	
	print("Selected "||hpf_qry("select count(distinct parent_sequence_key) from hpf_test.protein_seqs")||" proteins")
	return("hpf_test.protein_seqs")
}

domain_type_seqs = function(protein_table=NA, domain_type=NA, hpf=NA, confidence=0,class=NA,with_scop=F,two_filter=F){
	#print("ONE CALL")
	if(is.na(protein_table)){
		protein_table = protein_seqs(domain_type=domain_type)
	}
	query = "drop table if exists hpf_test.domain_seqs"
	#print(query)
	hpf_qry(query)
	
	domain_table="domain"
	query = "
	      select distinct d.domain_sequence_key as domain_sequence_key, d.domain_type, z.confidence
	      from "||domain_table||" d "
	if(!is.na(hpf)){
		query = query||"
			join bddb.filesystemOutfile o
			on o.sequence_key=d.domain_sequence_key
			and o.executable_key="||c(179,376)[hpf]
	}
  	query = query||"
	      join "||protein_table||" p 
	      	   on d.parent_sequence_key=p.parent_sequence_key" 
	if(with_scop){
		query = query||" join domain_sccs z on d.domain_sequence_key=z.domain_sequence_key "   
	}
	if(!is.na(domain_type)){
	   query = query||"
	      where 
	      	   d.domain_type in ('"||paste(domain_type,collapse="','")||"')
				and z.confidence >= "||confidence
	}
	if(two_filter){
		query=query||" and d.two_domain_filter='no'"
	}
   if(!is.na(class)){
	   classes = c()
	   for(c in class){
		   classes = c(classes,"sccs like '"||c||".%'")
	   }
	   query = query||" and ("||paste(classes,collapse=" or ")||")"
   }
   query = "create table hpf_test.domain_seqs "||query
   print(query)
   hpf_qry(query)
   
   query = "ALTER IGNORE TABLE hpf_test.domain_seqs ADD UNIQUE INDEX domain_sequence_key (domain_sequence_key)"
   #query = "create index domain_sequence_key on hpf_test.domain_seqs(domain_sequence_key)"
   #print(query) 
   hpf_qry(query)
   
   print("Selected "||hpf_qry("select count(distinct domain_sequence_key) from hpf_test.domain_seqs")||" domains")
   return("hpf_test.domain_seqs");

}

split = function(data){
      data_go = data[which(is.na(data$mf_acc)),]
      #data_i = data[which(!(data$mf_acc %in% NA) & (data$go_acc != data$mf_acc)),]
	  data_i = data[which((!(is.na(data$mf_acc)) & (is.na(data$go_acc))) | !(data$go_acc == data$mf_acc)),]
	  data_c = data[which(!(is.na(data$go_acc)) & (data$go_acc == data$mf_acc)),]
      data_go$labels = rep(1,length(data_go$mf_acc))
      data_i$labels = rep(-1,length(data_i$mf_acc))
      data_c$labels = rep(1,length(data_c$mf_acc))
      return(list(data_go,data_c,data_i))
}

bin_llr = function(data, breaks=c(0,-4,-8,-13), predictor="base_llr"){
	 bins = list()
	 names = c()
	 for(i in 1:(length(breaks)-1)){
	       bins = c(bins,list(data[which(data[[predictor]] <= breaks[i] & data[[predictor]] > breaks[i+1]),]))
	       names = c(names,""||breaks[i]||":"||breaks[i+1])
	 }
	 return(bins)
}

load = function(seq_table,base_min=-99,base_max=0,predictor_min=-99,bayes_table="functionPredictions.bayes_golite_062009_3",distance=0,unannotated=F,go_table="hddb_IEA.hddb_IEA_goLite_062009",random=F){
	
	# Setup temporary table
	hpf_qry("drop table if exists hpf_test.mfs")
	hpf_qry("drop table if exists allmf")
	hpf_qry("drop table if exists base")
	
	# Figure out all of the molecular functions and implied ancestors
	query = "create temporary table allmf
	select distinct t1.acc
	from "||seq_table||" as seq_table
	join hpf.domain d
		on seq_table.domain_sequence_key = d.domain_sequence_key
	join "||go_table||" f
		on d.parent_sequence_key = f.sequence_key
		and f.evidence_code in ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA','hpf_IEA' )
	join mygoLite_062009.term t2
		on f.acc = t2.acc
		and t2.term_type = 'molecular_function'
		and t2.is_obsolete = 0
	join mygoLite_062009.graph_path path
		on t2.id = path.term2_id
	join mygoLite_062009.term t1
		on path.term1_id=t1.id
		and t1.term_type='molecular_function'
		and t1.is_obsolete = 0
	"
	print(query)
	hpf_qry(query)
	hpf_qry("create index acc on allmf (acc)")
	
	# Get the base llr for all functions visible
	query = "create temporary table base
	select l.acc,l.metric
	from allmf m join functionTables.log_ratio_goLite_062009_3 l 
	on m.acc=l.acc and l.acc2 is NULL"
	hpf_qry(query)
	hpf_qry("create index acc on base (acc)")
	
	# Want to join with graph path to get all implied ancestor annotations
	# then join with the base_llr table
	query = "
	create table hpf_test.mfs ( 
	select DISTINCT 
	d.domain_sequence_key as domain_sequence_key,
	d.parent_sequence_key as parent_sequence_key,
	t1.id as mf_id, 
	t1.acc as mf_acc,
	t1.name as name,
	ifnull(l.metric,-12) as base_llr
	from "||seq_table||" as seq_table
	join hpf.domain d
		on seq_table.domain_sequence_key = d.domain_sequence_key
	join "||go_table||" f
		on d.parent_sequence_key = f.sequence_key
		and f.evidence_code in ( 'TAS','IDA','IMP','IGI','IPI','ISS','IEA','IC','IEP','NAS','ND','NR','RCA','hpf_IEA' )
	join mygoLite_062009.term t2
		on f.acc = t2.acc
		and t2.term_type = 'molecular_function'
		and t2.is_obsolete = 0
	join mygoLite_062009.graph_path path
		on t2.id = path.term2_id
	join mygoLite_062009.term t1
		on path.term1_id=t1.id
		and t1.term_type='molecular_function'
		and t1.is_obsolete = 0
	left outer join base l
		on t1.acc=l.acc
	)";
#base.metric as base_llr,
#join functionTables.log_ratio_goLite_062009 base
#on t2.acc = base.acc2 and base.acc is NULL

	print(sub("\t"," ",sub("\n"," ",query)))
	hpf_qry(query);

	query = "create index domain_sequence_key on hpf_test.mfs(domain_sequence_key)"
	#print(query)
	hpf_qry(query)
	query = "create index mf_acc on hpf_test.mfs(mf_acc)"
	#print(query)
	hpf_qry(query)

	
	#Crazy query joins GO annotations and Function Predictions into one table
	# With this we can see what GO has that we don't or the other way around.
	query = "
		select DISTINCT
			p.sequence_key as parent_sequence_key,
			ery.domain_sequence_key,
			seq_table.domain_type,
			seq_table.confidence as domain_confidence,
			ery.mf_acc,
			ery.name, 
			ery.p_llr, 
			ery.s_llr,
			ery.l_llr, 
			ery.ls_llr,
			ery.ps_llr,
			ery.pl_llr,
			ery.pls_llr, 
			ery.base_llr,
			ery.sccs,
			mfs.mf_id as mfs_mf_id,
			mfs.mf_acc as go_acc,
			mfs.name as go_name
		from "||seq_table||" as seq_table
		join "||bayes_table||" as ery
			on seq_table.domain_sequence_key = ery.domain_sequence_key
			and ery.base_llr <=  "||base_max||" 
			and ery.base_llr >=  "||base_min||"
			and ery.pls_llr >= "||predictor_min||"
		join hpf.protein p
			on ery.parent_sequence_key = p.sequence_key
		left outer join hpf_test.mfs as mfs 
			on ery.domain_sequence_key = mfs.domain_sequence_key 
			and ery.mf_acc = mfs.mf_acc"
 	
	if(unannotated){
		query = query||"	
		union distinct
		
		select DISTINCT
			p.sequence_key as parent_sequence_key,
			seq_table.domain_sequence_key,
			seq_table.domain_type,
			seq_table.confidence as domain_confidence,
			ery.mf_acc,
			ery.name,
			ery.p_llr, 
			ery.s_llr,
			ery.l_llr,
			ery.ls_llr,
			ery.ps_llr,
			ifnull(ery.pl_llr,-99),
			ifnull(ery.pls_llr,-99), 
			mfs.base_llr,
			ery.sccs,
			mfs.mf_id as mfs_mf_id,
			mfs.mf_acc as go_acc,
			mfs.name as go_name
		from "||seq_table||" as seq_table
		join hpf_test.mfs as mfs 
			on seq_table.domain_sequence_key=mfs.domain_sequence_key
			and mfs.base_llr <=  "||base_max||"
			and mfs.base_llr >=  "||base_min||"
		join hpf.protein p
			on mfs.parent_sequence_key=p.sequence_key
		left outer join "||bayes_table||" as ery
			on mfs.domain_sequence_key = ery.domain_sequence_key 
			and ery.mf_acc = mfs.mf_acc
			and ery.pls_llr >= "||predictor_min			
	}

	print(sub('\t',' ',sub('\n',' ',query)))
	#print(query)
	data = hpf_qry(query)

	if(random)
	{
		randquery = "select mf_acc from "||bayes_table||" where base_llr <= "||base_max||" and base_llr >= "||base_min||" order by RAND() limit "||length(data[[1]])
		print(randquery)
		rand_data = hpf_qry(randquery)
		print(rand_data)
		data$mf_acc = rand_data$mf_acc
	}

	return(data)

}

plotBayes <- function(seq_table,eval_table,functionFrom, predictors="s", annotations=FALSE, domains=FALSE, proteins=FALSE, break_points = c(-2,0,2,4,6,8,10,12,14,16,25,100), base_max=0, base_min=-99, bayes=TRUE)
{
	return(plotProbability(seq_table=seq_table, eval_table=eval_table, functionFrom=functionFrom, predictors=predictors, annotations=annotations, domains=domains,proteins=proteins, break_points=break_points, base_max=base_max, base_min=base_min,bayes=bayes));
}

#break_points = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.0)
#break_points = c(10**seq(0,1.05,.05)/10 -.1)
plotProbability <- function(seq_table, eval_table, functionFrom, predictors="s", annotations=FALSE, domains=FALSE, proteins=FALSE, break_points = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.0), base_max = 1.0, base_min = 0.0, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA', 'IC','IEP','NAS','ND','NR','RCA'), random=FALSE, rand_seq_table="rand_seq_map", offset=500, bayes=FALSE)
{
	plot_correct_vector = c();
	plot_correct_normalized_vector = c();
	plot_notcorrect_vector = c();
	plot_notcorrect_normalized_vector = c();

	total_domains = 0;
	total_proteins = 0;

	getFunction(code=code, functionFrom = functionFrom, create_table=TRUE, seq_table=seq_table);

	if(random)
	{
		randomize_seqs(seq_table,rand_seq_table=rand_seq_table, offset=offset);
	}
	query = ""
	for(i in seq(length(break_points)-1))
	{
		if(annotations)
		{
			query = paste("select count(*) from (");
		}
		else if(domains)
		{
			query = paste("select count(distinct domain_sequence_key) from (");
		}
		else if(proteins)
		{
			query = paste("select count(distinct parent_sequence_key) from (");
		}

		if(bayes)
		{
			query = paste(query, " select ery.domain_sequence_key, ery.mf_acc, ery.name, ery.p_llr, ery.s_llr, ery.ps_llr, ery.base_llr, mfs.mf_id as mfs_mf_id ");
		}
		else
		{
			query = paste(query, " select ery.domain_sequence_key, ery.mf_acc, ery.name, ery.p_prob, ery.s_prob, ery.ps_prob, ery.base_prob, mfs.mf_id as mfs_mf_id ");
		}

		if(proteins)
		{
			query = paste(query, ", d.parent_sequence_key");
		}
		query = paste(query, " from ",eval_table, " as ery ");

		if(random)
		{

			query = paste(query, " inner join ", rand_seq_table," as rand_seq_table on (rand_seq_table.seq_id = ery.domain_sequence_key) ");
			query = paste(query, " left outer join mfs on (rand_seq_table.rand_seq_id = mfs.seq_id and ery.mf_acc = mfs.mf_id) ");
			query = paste(query, " , ", seq_table," as seq_table");
		}
		else
		{
			query = paste(query, " left outer join mfs on (ery.domain_sequence_key = mfs.seq_id and ery.mf_acc = mfs.mf_acc) ");
			query = paste(query, " , ", seq_table," as seq_table");
		}

		if(proteins)
		{
			query = paste(query, ", hpf.domain as d", sep="");
		}

		query = paste(query, " where seq_table.seq_id = ery.domain_sequence_key ");
		#query = paste(query, " where ");
		if(proteins)
		{
			query = paste(query, " and d.domain_sequence_key = ery.domain_sequence_key");
		}

		if(bayes)
		{
			query = paste(query, " and ery.base_llr <= ",base_max, " and ery.base_llr >= ",base_min);
			query = paste(query, " and ery.",predictors,"_llr >= ", break_points[i], " and ery.",predictors,"_llr < ",break_points[i+1], sep="");
		}
		else
		{
			query = paste(query, " and ery.base_prob <= ",base_max, " and ery.base_prob >= ",base_min);
			query = paste(query, " and ery.",predictors,"_prob >= ", break_points[i], " and ery.",predictors,"_prob < ",break_points[i+1], sep="");
		}
		query = paste(query, ") as tmp");

		correct_query = paste(query," where mfs_mf_id is not NULL");
		print(correct_query);
		correct = hpf_qry(correct_query)
		print(paste("correct:",correct));

		if(annotations)
		{
			notcorrect_query = paste(query," where mfs_mf_id is NULL");
			print(notcorrect_query);
			notCorrect = hpf_qry(notcorrect_query)
		}
		else if(proteins)
		{
			print(query);
			total_proteins = hpf_qry(query);
			notCorrect = total_proteins - correct;
		}
		else if(domains)
		{
			print(query);
			total_domains = hpf_qry(query);
			notCorrect = total_domains - correct;
		}

		print(paste("notCorrect:",notCorrect));

		plot_correct_vector <- c(plot_correct_vector, correct[1,1]);
		plot_notcorrect_vector <- c(plot_notcorrect_vector, notCorrect[1,1]);
		plot_correct_normalized_vector <- c(plot_correct_normalized_vector, correct[1,1]/(correct[1,1]+notCorrect[1,1]));
		plot_notcorrect_normalized_vector <- c(plot_notcorrect_normalized_vector, notCorrect[1,1]/(correct[1,1]+notCorrect[1,1]));
	}

		plot_matrix = rbind(plot_correct_vector, plot_notcorrect_vector);
		plot_normalized_matrix = rbind(plot_correct_normalized_vector, plot_notcorrect_normalized_vector);

		print(plot_matrix);
		print(plot_normalized_matrix);
		#plotProbability_barplot(plot_matrix, break_points);

		return(plot_matrix);
}

plot_specificity <- function(eval_table, predictors="s", background_break_points = c(-10, -8, -6,-4, -2, 0), pred_prob_min = 0)
{
	#create temporary table scop_01_001 select distinct seq_id  from evaluation_results_scop where ps_prob > .5 and base_prob <.01 and base_prob >=.001 and seq_id not in (select seq_id from scop_001_0001);

	for(i in seq(length(background_break_points)-1))
	{
		drop_qry = paste("drop table if exists seqs_",i,"_",i+1,sep="");
		print(drop_qry);
		hpf_qry(drop_qry);

		create_qry = paste("create temporary table seqs_",i,"_",i+1,sep="");
		create_qry = paste(create_qry, " select distinct domain_sequence_key from ",eval_table);
		create_qry = paste(create_qry, " where ",predictors,"_llr >",pred_prob_min," and base_llr <",sprintf("%f",background_break_points[i+1])," and base_llr >=",sprintf("%f",background_break_points[i]), sep="");
		#print(background_break_points[i])
		#print(background_break_points[i+1])
		for(j in seq(length(background_break_points)-1))
		{
			if(j < i)
			{
				create_qry = paste(create_qry," and domain_sequence_key not in (select domain_sequence_key from seqs_",j,"_",j+1,")",sep="");
			}
		}

		print(create_qry);
		hpf_qry(create_qry);
	}
	#select count(distinct seq_id) from scop_001_0001
	spcfy <- c();
	spcfy_names <- c();
	for(i in seq(length(background_break_points)-1))
	{
		count_qry = paste("select count(distinct domain_sequence_key) from seqs_",i,"_",i+1,sep="");
		print(count_qry);
		cnt = hpf_qry(count_qry);
		spcfy <- c(spcfy,cnt[,1]);
		spcfy_names <- c(spcfy_names, paste(background_break_points[i+1],"_",background_break_points[i],": ", cnt[,1], sep=""));
	}
	names(spcfy) = spcfy_names;
	pie(spcfy);
}

create_benchmark_sequences <- function(benchmark_db, sf_data_table="mcm_data", table_name, MF=TRUE, BP=FALSE, CC=FALSE, code=c('TAS','IDA','IMP','IGI','IPI','ISS','IEA', 'IC','IEP','NAS','ND','NR','RCA'))
{
	go_evidence = paste('\'',code,'\'', collapse=',', sep='');

	#create table arabidopsis_fr_mf_seqs select distinct fd.domain_sequence_key as seq_id from fr_data as fd, function as f, mygo.term as t where fd.parent_sequence_key = f.sequence_key and f.acc = t.acc and t.term_type = "molecular_function";
	#select distinct fd.*, t.id, t.name, f.code, t2.id, t2.name, f2.code from fr_data as fd, function as f, mygo.term as t, function as f2, mygo.term as t2 where fd.parent_sequence_key = f.sequence_key and fd.parent_sequence_key = f2.sequence_key and f.acc = t.acc and f2.acc = t2.acc and t.term_type = "molecular_function" and t2.term_type = "biological_process"

	sql_qry = paste("create table ", benchmark_db,".",table_name, sep="");
	sql_qry = paste(sql_qry, " select distinct d.domain_sequence_key as seq_id ");
	sql_qry = paste(sql_qry, " from ", benchmark_db,".",sf_data_table," as d", sep="");
	
	if(MF)
	{
		sql_qry = paste(sql_qry, ", " , benchmark_db , ".function as f", sep="");
		sql_qry = paste(sql_qry, ", mygo.term as t");
	}
	if(BP)
	{
		sql_qry = paste(sql_qry, ", ",benchmark_db,".function as f2", sep="");
		sql_qry = paste(sql_qry, ", mygo.term as t2");
	}
	if(CC)
	{
		sql_qry = paste(sql_qry, ", ",benchmark_db,".function as f3", sep="");
		sql_qry = paste(sql_qry, ", mygo.term as t3");
	}

	sql_qry = paste(sql_qry, " where ");

	and_flag = FALSE;
	if(MF)
	{
		sql_qry = paste(sql_qry, " d.parent_sequence_key = f.sequence_key and f.acc = t.acc and t.term_type = 'molecular_function' and f.code in (", go_evidence, ") and t.acc <> 'GO:0003674'");
		and_flag = TRUE;
	}
	if(BP)
	{
		if(and_flag)
		{ sql_qry = paste(sql_qry, " and ");}
		sql_qry = paste(sql_qry, " d.parent_sequence_key = f2.sequence_key and f2.acc = t2.acc and t2.term_type = 'biological_process' and f.code in (", go_evidence, ") and t2.acc <> 'GO:0008150'");
		and_flag = TRUE;
	}
	if(CC)
	{
		if(and_flag)
		{ sql_qry = paste(sql_qry, " and ");}
		sql_qry = paste(sql_qry, " d.parent_sequence_key = f3.sequence_key and f3.acc = t3.acc and t3.term_type = 'cellular_component' and f.code in (", go_evidence, ") and t3.acc <> 'GO:0005575'");
	}

	print(sql_qry);
	hpf_qry(sql_qry);
	
	alter_qry = paste("alter table ", benchmark_db,".",table_name, " add index(seq_id)", sep="");
	print(alter_qry);
	hpf_qry(alter_qry);
		
}

normalize_view <- function(plot_matrix)
{
	correct = plot_matrix[1,]/(plot_matrix[1,]+plot_matrix[2,]);
	not_correct = plot_matrix[2,]/(plot_matrix[1,]+plot_matrix[2,]);

	print(plot_matrix);
	print(rbind(correct, not_correct));
}

plotProbability_lineplot <- function(plot_matrix, break_points, redraw=FALSE, color=1, bayes=FALSE)
{
	plot_normalized_matrix = plot_matrix[1,]/(plot_matrix[1,]+plot_matrix[2,]);

	if(redraw)
	{
		if(bayes)
		{
			plot(0:1.05, xlab = "log likelihood ratio", ylab="percent of total annotations in bin", type="n", asp=20);
		}
		else
		{
			plot(0:1,0:1.05, xlab = "probability", ylab="percent of total annotations in bin", type="n");
		}
	}

	lines(break_points, plot_normalized_matrix, type="b", col=color, cex=.5);
	for(i in seq(length(break_points)))
	{
		if((plot_matrix[1,i]+plot_matrix[2,i]) < 3000)
		{
			points(break_points[i], plot_normalized_matrix[i], type="b", cex=((plot_matrix[1,i]+plot_matrix[2,i])/50), col=color);
		}
	}

	#kdrew: start at 2.5 and add 3 
	x_pos_correct = 1.5;
	#x_pos_not_correct = 2.5;
	for(i in seq(length(break_points)-1))
	{
		text(x_pos_correct, plot_normalized_matrix[i]+.01, plot_matrix[1,i]);
		#text(x_pos_not_correct, plot_normalized_matrix[2,i]+.01, plot_matrix[2,i]);
		x_pos_correct = x_pos_correct + 3;
		#x_pos_not_correct = x_pos_not_correct + 3;
	}
}
plotProbability_barplot <- function(plot_matrix, break_points)
{
	names_vec = c();
	for(i in seq(length(break_points)-1))
	{
		name = paste(format(break_points[i], digits=2),"-\n",format(break_points[i+1],digits=2));
		names_vec = c(names_vec,name);
	}
	plot_normalized_matrix = plot_matrix[1,]/(plot_matrix[1,]+plot_matrix[2,]);
	plot_normalized_matrix = rbind(plot_normalized_matrix, (plot_matrix[2,]/(plot_matrix[1,]+plot_matrix[2,])));
	#barplot(plot_normalized_matrix, beside=TRUE, col=c("green","red"), names.arg=format(break_points[-1], digits=1), xlab = "probability", ylab="percent of total annotations in bin", ylim=c(0,1.05));
	#barplot(plot_normalized_matrix, beside=TRUE, col=c("steelblue","gray60"), names.arg=format(break_points[-1], digits=1), xlab = "probability (log binned)", ylab="percent of total annotations in bin", ylim=c(0,1.05));
	barplot(plot_normalized_matrix, beside=TRUE, col=c("steelblue","gray60"), names.arg=names_vec, xlab = "probability (log binned)", ylab="percent of total annotations in bin", ylim=c(0,1.05));
	title(sub="blue = correct, gray = incorrect (or not annotated yet)")

	#kdrew: start at 2.5 and add 3 
	x_pos_correct = 1.5;
	x_pos_not_correct = 2.5;
	for(i in seq(length(break_points)-1))
	{
		text(x_pos_correct, plot_normalized_matrix[1,i]+.01, plot_matrix[1,i]);
		text(x_pos_not_correct, plot_normalized_matrix[2,i]+.01, plot_matrix[2,i]);
		x_pos_correct = x_pos_correct + 3;
		x_pos_not_correct = x_pos_not_correct + 3;
	}
}


go_distribution <- function(seq_ids, eval_table, functionFrom, background=TRUE)
{
	#select t2.name, t2.acc, count(distinct sequence_key) as cnt, count(distinct sequence_key)/2597 as percent_background from function as f, mygo.term as t, mygo.graph_path as gp, mygo.term as t2 where t.acc = f.acc and gp.term2_id = t.id and gp.term1_id = t2.id and t2.term_type = "molecular_function" group by t2.id order by cnt DESC
	query = "";
	if(background)
	{
		total_query = paste("select @total_seqs := count(distinct sequence_key) from ",functionFrom,".function as f, mygo.term as t where t.acc = f.acc and t.term_type = 'molecular_function'", sep="");
		query = paste(query, "select t2.name, t2.acc, count(distinct sequence_key) as cnt, count(distinct sequence_key)/@total_seqs as percent_background ");
		query = paste(query, " from ",functionFrom,".function as f, mygo.term as t, mygo.graph_path as gp, mygo.term as t2 ", sep="");
		query = paste(query, " where t.acc = f.acc and gp.term2_id = t.id and gp.term1_id = t2.id and t2.term_type = 'molecular_function' ");
		query = paste(query, " group by t2.id order by cnt DESC");

		print(total_query);
		print(query);
		hpf_qry(total_query);
		result = hpf_qry(query);
	}

	return(result);
}

plotPredictions <- function(seq_ids, lines=FALSE, eval_table = "evaluation_results_PDB", functionFrom="scopfold", rank_thres = 100, predictors="s")
{
	flag_new = TRUE;
	lapply(seq_ids,function(seq_id)
	{
		query = paste("select mf_id, ",predictors,"_prob as predictors_prob, base_prob, ",predictors,"_rank as predictors_rank, base_rank,  ",predictors,"_rank_diff_rank as predictors_rank_diff_rank, ",predictors,"_prob_diff_rank as predictors_prob_diff_rank,if(",predictors,"_prob = base_prob,0,((1.0 - base_prob)- (1.0 - ",predictors,"_prob))/(1.0-base_prob)) as prob_score, (base_rank-",predictors,"_rank)/base_rank as rank_score from ",eval_table," where seq_id = ",seq_id," and ",predictors,"_rank <= ",rank_thres,sep="")

		#kdrew: mfs is a table created by getFunctionScopfold
		query2 = paste("select mf_id, ",predictors,"_prob as predictors_prob, base_prob, ",predictors,"_rank as preditors_rank, base_rank,  ",predictors,"_rank_diff_rank as predictors_rank_diff_rank, ",predictors,"_prob_diff_rank as predictors_prob_diff_rank,if(",predictors,"_prob = base_prob,0,((1.0 - base_prob)- (1.0 - ",predictors,"_prob))/(1.0-base_prob)) as prob_score, (base_rank-",predictors,"_rank)/base_rank as rank_score from ",eval_table, " as erp where seq_id = ",seq_id," and ",predictors,"_rank <= ",rank_thres," and mf_id in (select mf_id from mfs)", sep="")

		#kdrew: mfs is a table created by getFunctionScopfold
		query3 = paste("select mf_id, ",predictors,"_prob as predictors_prob, base_prob, ",predictors,"_rank as predictors_rank, base_rank,  ",predictors,"_rank_diff_rank as predictors_rank_diff_rank, ",predictors,"_prob_diff_rank as predictors_prob_diff_rank,if(",predictors,"_prob = base_prob,0,((1.0 - base_prob)- (1.0 - ",predictors,"_prob))/(1.0-base_prob)) as prob_score, (base_rank-",predictors,"_rank)/base_rank as rank_score from ", eval_table, " as erp where seq_id = ",seq_id," and ",predictors,"_rank <= ",rank_thres," and mf_id in (select mf_id from mfs)", sep="")

		print(query);
		allData = hpf_qry(query)

		query_path = paste("select * from mygo.graph_path as gp, mfs where mfs.mf_id = gp.term2_id and distance = 1 and term1_id <> 1");

		#getFunctionScopfold(seq_id, create_table=TRUE);
		getFunction(seq_id, functionFrom = functionFrom, create_table=TRUE);
		print(query2);
		trueData_high = hpf_qry(query2)
		paths = hpf_qry(query_path);

		#getFunctionScopfold(seq_id, code=c('IC','IEP','NAS','ND','NR','RCA'), create_table=TRUE);
		getFunction(seq_id, code=c('IC','IEP','NAS','ND','NR','RCA'), functionFrom=functionFrom, create_table=TRUE);
		print(query3);
		trueData_low = hpf_qry(query3)
		paths = rbind(paths,hpf_qry(query_path));


		if(flag_new)
		{
			plot(allData[,"rank_score"],allData[,"prob_score"],pch=".", xlim=c(-1,1),ylim=c(-1,1))
			#plot(allData[,"s_prob_diff_rank"],allData[,"s_rank_diff_rank"],pch=".")
			#plot(allData[,"s_rank"],allData[,"s_prob"],pch=".")
			flag_new <<- FALSE;
		}
		else
		{
			points(allData[,"rank_score"],allData[,"prob_score"], pch=".")
			#points(allData[,"s_prob_diff_rank"],allData[,"s_rank_diff_rank"], pch=".")
			#points(allData[,"s_rank"],allData[,"s_prob"], pch=".")
		}
		if(length(trueData_low) > 0)
		{
			points(trueData_low[,"rank_score"],trueData_low[,"prob_score"], pch=20, col="blue")
			#points(trueData_low[,"s_prob_diff_rank"],trueData_low[,"s_rank_diff_rank"], pch=20, col="blue")
			#points(trueData_low[,"s_rank"],trueData_low[,"s_prob"], pch=20, col="blue")
		}


		if(length(trueData_high) > 0)
		{
			points(trueData_high[,"rank_score"],trueData_high[,"prob_score"], pch=20, col="red")
			#points(trueData_high[,"s_prob_diff_rank"],trueData_high[,"s_rank_diff_rank"], pch=20, col="red")
			#points(trueData_high[,"s_rank"],trueData_high[,"s_prob"], pch=20, col="red")
		}


		if(lines)
		{
			apply(paths, 1, function(path)
			{
				#print(path[2]);
				scores_vec = c(allData[allData$mf_id == path["term1_id"],"rank_score"],allData[allData$mf_id == path["term2_id"],"rank_score"]);
				#scores_vec = c(allData[allData$mf_id == path["term1_id"],"s_prob_diff_rank"],allData[allData$mf_id == path["term2_id"],"s_rank_diff_rank"]);
				#scores_vec = c(allData[allData$mf_id == path["term1_id"],"s_rank"],allData[allData$mf_id == path["term2_id"],"s_prob"]);
				s_prob_vec = c(allData[allData$mf_id == path["term1_id"],"prob_score"],allData[allData$mf_id == path["term2_id"],"prob_score"]);
				#s_prob_vec = c(allData[allData$mf_id == path["term1_id"],"s_prob_diff_rank"],allData[allData$mf_id == path["term2_id"],"s_rank_diff_rank"]);
				#s_prob_vec = c(allData[allData$mf_id == path["term1_id"],"s_rank"],allData[allData$mf_id == path["term2_id"],"s_prob"]);

				#print(scores_vec);
				#print(s_prob_vec);
				lines(scores_vec,s_prob_vec);
			
			});
		}

		points(allData[allData$mf_id == 1181,"rank_score"],allData[allData$mf_id == 1181,"prob_score"], pch=19, col="green")
		#points(allData[allData$mf_id == 1181,"s_prob_diff_rank"],allData[allData$mf_id == 1181,"s_rank_diff_rank"], pch=19, col="green")
		#points(allData[allData$mf_id == 1181,"s_rank"],allData[allData$mf_id == 1181,"s_prob"], pch=19, col="green")
	});
}

