library(RMySQL)
# Just hold a connection for simplicity
if(!exists("hpf")){
	hpf <- dbConnect(MySQL(), dbname="hpf", password="patrick_nyu")
}

# Creating a string concatenation operator
"||" <- function(...) UseMethod("||")
"||.default" <- .Primitive("||")
"||.character" <- function(...) paste(...,sep="")

# Example: "abc" || "def" || "ghi" 

split = function(data, since_solved=T, adjusted=F){
	data_c = correct(data, since_solved=since_solved, adjusted=adjusted)
	data_c$labels = rep(1,length(data_c$id))
	data_i = data[which(!(data$id %in% data_c$id)), ]
	data_i$labels = rep(-1,length(data_i$id))
	return(list(data_c,data_i))
}

#' Grab a subset of the data with the given ids
#' @param data the dataframe
#' @param ids the vector of indices
#' @returnType data.frame
#' @return the data subset
#' @author patrick
#' @export
id_subset = function(data, ids){
	sub = which(data$sequence_key %in% ids)
	data_sub = data[sub,]
	return(data_sub)
}

subset_number = function(top, number=5, since_solved=FALSE){
	top_indices = c()
	if(since_solved){
		the_set = unique(top$cluster_id)
	} else {
		the_set = unique(top$sequence_key)
	}
	for(i in the_set){
		if(since_solved){
			which_ones = which(top$cluster_id %in% i)
		} else {
			which_ones = which(top$sequence_key %in% i)
		}
			
		#print("Length "||length(unique(which_ones)))
		#If there arent "number" make sure you do a min or it will duplicate values
		top_indices=c(top_indices,which_ones[1:min(number, length(which_ones))])
	}
	
	return(top[top_indices,])

}

predict_alpha = function(data){
	return(data[which(data$prediction_percent_beta <= 0.15 & data$prediction_percent_alpha >= 0.15),])	
}

predict_beta = function(data){
	return(data[which(data$prediction_percent_beta >= 0.15 & data$prediction_percent_alpha <= 0.15),])	
}

predict_other = function(data){
	return(data[which((data$prediction_percent_beta >= 0.15 & data$prediction_percent_alpha >= 0.15) | (data$prediction_percent_beta <= 0.15 & data$prediction_percent_alpha <= 0.15)),])	
}

#' Substrings and removes erroneous sccs information.  If you pass a data.frame, replaces the "sccs" column with all names formatted.
#' @param character or data.frame
#' @returnType character or data.frame
#' @return The 3 part superfamily id.
#' @author patrick
#' @export
superfamily = function(object){
	if(is.data.frame(object)){
		experiment_sccs = c()
		for(name in object$experiment_sccs){
			experiment_sccs = c(experiment_sccs, superfamily(name))
		}
		object$experiment_sccs = experiment_sccs
		return(object)
	} else if (is.character(object)) {
		return(paste(unlist(strsplit(object,"\\."))[1:3], collapse="."))
	}
}

rmsd = function(hpf_id=2, ids=NA){
	in_id = ""
	if(!is.na(ids)){
		in_id = " and sequence_key in ("||paste(ids, collapse=", ")||")"
	}
	query = "select r.* from hpf2.rmsd r where r.analysis_key="||hpf_id||""||in_id||""
	print(query)
	rmsd = dbGetQuery(hpf,query)
}

#' Searches for true sccs for domains in data, filtering out scores on domains that don't have this class.  The default blank string for class should give back all domains for which we know the correct sccs.
#' @param data the dataframe
#' @param class the class subset you want
#' @returnType data.frame
#' @return data for domains whose real classification is in this class
#' @author patrick
#' @export
scop_class = function(data, class="", since_solved=F){
	indices=c()
	for(domain in unique(data$sequence_key)){
		if(!is.na(domain)){
			if(since_solved){
				query = "select concat(s.sccs,\".1\") from since_solved.since_solved s join since_solved.cdhit_clstr c on s.cluster_id=c.cluster_id where c.identifier="||domain
			} else {
				query = "select a.sccs from hpf2.benchmark b join pdb.astral a on b.astral_key=a.id where b.sequence_key="||domain||";" 
			}
			correct_sccs = dbGetQuery(hpf, query)
			# Check true sccs with the regular expression
			#print(correct_sccs)
			if(length(grep(class||"\\..*\\..*\\..*", correct_sccs)) > 0){
				indices = c(indices, which(data$sequence_key %in% domain))
			} else {
				#print("no")
			}
		}
	}
	return(data[indices,])
}

#' 
#' @param filter_id 
#' @returnType data.frame
#' @return the top5 cluster index scores per domain sequence id for the given filter.
#' @author patrick
#' @export
top = function(filter_id=2, experiment=NA, number=5, hpf_id=2, ids=NA, adjusted=FALSE, scale=0, since_solved=FALSE,date="2005-01-01"){

    if(since_solved){
	  filter_id=4
	  query = since_solved_query(adjusted=adjusted, scale=scale, filter_id=filter_id,date=date)
	} else {
	  query = experiment_query(filter_id=filter_id, experiment=experiment, hpf_id=hpf_id, adjusted=adjusted, scale=scale)
	}
	print(query)
	top = dbGetQuery(hpf,query)
	
	if(!is.na(ids)){
		top = subset(top, ids)
	}
	top = superfamily(top)
	top = subset_number(top, number=number, since_solved=since_solved)
	return(top)
}

experiment_query = function(filter_id=2, experiment=NA, hpf_id=2, adjusted=FALSE, scale=0){
	where = "where e.id="||experiment
	if(identical(hpf_id,2)){
		if(adjusted){
			table = "hpf2.mcmData_adjusted"
			if(!is.na(scale)){
				where=where||" and m.scale="||scale
			}
		} else {
			table = "hpf2.mcmData"
		}
		where = where||" and m.filter_id="||filter_id
	} else if (identical(hpf_id,1)){
		table = "hpf.mcmData"
	}	
	query="select m.*,length(s.sequence) as length from experiment e join protein p join domain d join "||table||" m join hpf.sequence s on e.id=p.experiment_key and p.sequence_key=d.parent_sequence_key and d.domain_sequence_key=m.sequence_key and m.sequence_key=s.id "||where||" order by sequence_key, probability desc;"
	return(query);
	
}

since_solved_query = function(filter_id=4, adjusted=FALSE,scale=0,date='2005-01-01'){
	if(adjusted){
	  table="hpf2.mcmdata_adjusted"
	  scale_q=" and scale="||scale
	  query = "SELECT m.*,s.cluster_id FROM  since_solved.since_solved s join "||table||" m on s.cluster_id = m.sequence_key where filter_id="||filter_id||scale_q||" and ascession_date>DATE('"||date||"') order by s.cluster_id, probability desc;"
	  
        } else {
	  table="hpf.mcmData"
	  scale_q=""
	  query = "SELECT m.*,s.cluster_id,o.* FROM  since_solved.since_solved s join since_solved.cdhit_clstr c on s.cluster_id=c.cluster_id join "||table||" m on c.identifier = m.sequence_key left outer join since_solved.mammoth o on m.structure_key=o.prediction where ascession_date>DATE('"||date||"') order by s.cluster_id, probability desc;"
	  
	}

	
	return(query)
}

label = function(ss, adjusted=F,since_solved=T){
	ss_c = mcm_correct(ss,adjusted=adjusted,since_solved=since_solved)
	ss_c$labels = rep(1,length(ss_c$id))
	ss_i = ss[which(!(ss$id %in% ss_c$id)), ]
	ss_i$labels = rep(-1,length(ss_i$id))
	return(rbind(ss_c,ss_i))
}


#' Generate a subset of dataframe with only correct scores.  Correct SCCS pulled from pdb benchmark and since_solved data tables;
#' @param data 
#' @returnType 
#' @return 
#' @author patrick
#' @export
mcm_correct = function(data, adjusted=FALSE, since_solved=FALSE){
	correct_indices = c()
	for(domain in unique(data$sequence_key)){
		if(!is.na(domain)){
			if(since_solved){
				query = "select s.sccs from since_solved.since_solved s join since_solved.cdhit_clstr cl on s.cluster_id=cl.cluster_id where cl.foreign_key="||domain
				#}
				correct_sccs = dbGetQuery(hpf,query)
				#str(correct_sccs)
			} else {
				query = "select a.sccs from hpf2.benchmark b join pdb.astral a on b.astral_key=a.id where b.sequence_key="||domain
				pdb_sccs = dbGetQuery(hpf, query)
				query = "select s.sccs from since_solved.since_solved s join since_solved.cdhit_clstr cl on s.cluster_id=cl.cluster_id where cl.foreign_key="||domain
				solved_sccs = dbGetQuery(hpf,query)
				correct_sccs = rbind(pdb_sccs,solved_sccs)
			}

			correct_sccs = superfamily(correct_sccs$sccs)
			indices = unique(which(data$sequence_key %in% domain & data$experiment_sccs %in% correct_sccs))
			#print(data[which(data$sequence_key %in% domain), ]$experiment_sccs)
			
			#print(paste(domain, correct_sccs))
			#print("number correct: "||length(indices))
			#print("indices: "||indices)
			if(length(indices)){
				#print(data$sccs[,indices])
				correct_indices = c(correct_indices, indices)
			}
		}		
	}
	return(data[correct_indices,])
}

correct_rmsd = function(data, rmsd){
	cor = correct(data)
	all_indices = c()
	for(domain in unique(cor$sequence_key)){
		good_clusters = cor[which(cor$sequence_key %in% domain),]$cluster_center_index
		
		indices = which(rmsd$sequence_key %in% domain & rmsd$cluster_center_index %in% good_clusters)
		all_indices = c(all_indices,indices)
	}
	return(rmsd[all_indices,])
}

#top_hpf1[["mcmdata_key"]]
average = function(data, index, column){
	ident = c()
	m = c()
	
	sequences = unique(data[[index]])
	for(i in sequences){
		indices = which(data[[index]] %in% i)
		ident = c(ident, i)
		print("sequence "||i||" mean "||mean(data[indices,][[column]])||" data "||data[indices,][[column]])
		m = c(m, mean(data[indices,][[column]]))
	}
	return(list(index=ident, column=m))
}
	
common = function(left, right){
	left_indices = c()
	right_indices = c()
	
	for(i in (1:length(left$id))){
		domain = left[i,]$sequence_key
		index = left[i,]$cluster_center
		
		l = i
		r = which(right$sequence_key %in% domain & right$cluster_center_index %in% index)
		
		if(length(r)>=1){
			right_indices = c(right_indices, r)
			left_indices = c(left_indices, l)
		} else {
			print(r)
		}
	}
	return(list(left=left[unique(left_indices),], right=right[unique(right_indices),]))
}

scatter_class = function(data, since_solved=F, adjusted=F,	classes = c("a","b","c","d"), colors = c("red","blue","green","orange", "black")){
	data_correct = correct(data, since_solved=since_solved, adjusted=adjusted)
	print("Correct "||length(unique(data_correct$sequence_key)))
	ids = c()
	for(i in (0:length(classes)+1)){
		print(""||i||":"||classes[i]||" "||colors[i])
		if(i <= length(classes)){
			data_class = class(data, class=classes[i], since_solved=since_solved)
			print(""||classes[i]||" "||length(unique(data_class$sequence_key)))
		} else {
			data_class = data[which(!data$id %in% ids), ]
		}
		ids = c(ids,data_class$id)
		if(i != 1){
			par(new=T)
		}
		data_class_correct = data_class[which(data_class$id %in% data_correct$id), ]
		print("Correct "||classes[i]||" "||length(unique(data_class_correct$sequence_key)))
		data_class_incorrect = data_class[which(!data_class$id %in% data_class_correct$id), ]
		print("Incorrect "||classes[i]||" "||length(unique(data_class_incorrect$sequence_key)))
		
		plot(data_class_incorrect$prediction_percent_alpha, data_class_incorrect$prediction_percent_beta, xlim=c(0,1), ylim=c(0,1), col=colors[i])
		par(new=T)
		plot(data_class_correct$prediction_percent_alpha, data_class_correct$prediction_percent_beta, xlim=c(0,1), ylim=c(0,1), col=colors[i], pch=23)
		
		#pch=23
	}
	par(new=T)
	plot(c(0,1),c(.15,.15), type="l", xlim=c(0,1),ylim=c(0,1), col="red")
	par(new=T)
	plot(c(.15,.15),c(0,1), type="l", xlim=c(0,1),ylim=c(0,1), col="blue")
	
}

scatter_rmsd = function(data, rmsd, main="MCM/RMSD Scatter Plot"){

	
	common = function(data, rmsd){
		rms = c()
		prob = c()
		domains = c()
		clusters = c()
		for(i in (1:length(data$id))){
			domain = data[i,]$sequence_key
			index = data[i,]$cluster_center_index
			
			probability = data[i,]$probability
			if(length(probability)>1){
				print("PLonger "||length(probability))
			}
			
			rmsd_score = rmsd[which(rmsd$sequence_key %in% domain & rmsd$cluster_center_index %in% index),]$rmsd
			if(length(rmsd_score) > 1){
				print("Longer "||length(rmsd_score))
				print(unique(rmsd_score))
			}
			
			
			if(length(rmsd_score)==1 & length(probability)==1){
				rms = c(rms, rmsd_score)
				prob = c(prob, probability)
				domains = c(domains, domain)
				clusters = c(clusters, index)
				#print(""||rms||" "||prob)
			}
		}
		return(list(rms=rms, prob=prob, sequence_key=domains, cluster_center_index=clusters))
	}
	
	results = common(data, rmsd)
	rms = results$rms
	prob = results$prob
	print(""||length(rms)||" points")
	plot(y=prob, x=rms, ylim=rev(range(prob)),xlab="RMSD", ylab="MCM Score", main=main||": "||length(unique(results$sequence_key))||" sequences", col="blue")
	all = list(mcm=prob,rmsd=rms, sequence_key=results$sequence_key, cluster_center_index=results$cluster_center_index)
	
	c_data = correct(data)
	c_rmsd = correct_rmsd(data, rmsd)
	results = common(c_data, c_rmsd)
	rms = results$rms
	prob = results$prob
	print(""||length(rms)||" points correct")
	points(y=prob,x=rms, col="seagreen" ,pch=19)
	cor = list(mcm=prob, rmsd=rms, sequence_key=results$sequence_key, cluster_center_index=results$cluster_center_index)
	
	return(list(all=all, correct = cor))

}

#' Run correct() on the data and plot a histogram overlaying correct scores with hist_overlay().
#' @param data 
#' @returnType NA
#' @return NA
#' @author patrick
#' @export
correct_hist = function(l, xlab="MCM Score", ylab="Count", title="Scop Benchmark Scores", col=colors()[c(377,401,516,574)], breaks=20, since_solved=FALSE, adjusted=FALSE){
	overlay_list = list()  
	for(i in (1:length(l))){
		data = l[[i]]
		good = correct(data, since_solved=since_solved, adjusted=adjusted)
		overlay_list = c(overlay_list, list(data,good))
	}
  	mcm_overlay(overlay_list, xlab=xlab, ylab=ylab, title=title, labels=TRUE, col=col, breaks=breaks)
}

compare = function(five, top, col=colors()[c(377,401,516,574)], main="HPF1/HPF2 Scop Benchmark Comparison"){
	correct_mhist(five, col=col[1:2], main=main)
	correct_mhist(top, col=col[3:4], add=TRUE, axes=FALSE, labels=FALSE, main=NA)
}

correct_mhist = function(l, main="Scop Benchmark", col=colors()[c(377,401)], add=FALSE, axes=TRUE, labels=TRUE, xlab="MCM Score", ylab="Count",freq=TRUE){
	h = hist(l[[1]]$probability, breaks=20, plot=FALSE)
	num_domains = length(unique(l[[1]]$sequence_key))
	if(!is.na(main)){
		main=main||": "||num_domains||" domains"
	}
	
	# Plot yields
	l_probability = list()
	for(i in (1:length(l))){
		l_probability = c(l_probability, list(l[[i]]$probability))
	}
	multhist(l_probability, main=main, beside=TRUE, col=col[1], breaks=h$breaks, add=add, axes=axes, xlab=xlab, ylab=ylab, labels=labels,freq=freq)

	# Now calculate correct
	l_correct = list()
	for(i in (1:length(l))){
		l_correct = c(l_correct, list(correct(l[[i]])$probability))
	}
	multhist(l_correct, beside=TRUE, col=col[2], breaks=h$breaks, add=TRUE, axes=FALSE, labels=FALSE,freq=freq)
	
}

correct_rmsd_hist = function(rmsd, data, xlab="RMSD", ylab="Count", title="RMSD to Native Structure"){
	good = correct_rmsd(data, rmsd)
	rmsd_overlay(rmsd$rmsd, good$rmsd, xlab=xlab, ylab=ylab, title=title||" "||length(unique(rmsd$sequence_key))||" Domains")
}

correct_rmsd_endpsi_hist = function(rmsd, data, xlab="Percent Sequence Alignment", ylab="Count", title="Mammoth Best Substructure Sequence Alignment"){
	good = correct_rmsd(data, rmsd)
	rmsd_overlay(rmsd$end_psi, good$end_psi, xlab=xlab, ylab=ylab, title=title||" "||length(unique(rmsd$sequence_key))||" Domains")
}

#' Overlays two histograms.
#' multhist(l) is an alternative side-by-side
#' @param best 
#' @param good 
#' @returnType 
#' @return 
#' @author patrick
#' @export
mcm_overlay = function(l, labels=FALSE, xlab="MCM Score", ylab="Count", title="", col=NA, breaks=20){
	if(is.na(col)){
		col=sample(colors(), length(l))
	}
	
	i=1
	hists = list()
	h = hist(l[[1]]$probability, xlab=xlab, ylab=ylab, main=title||" "||length(unique(l[[1]]$sequence_key))||" Domains", breaks=breaks, col=col[i])
	hists = c(hists,list(h))
	while(i<length(l)){
		i=i+1
		h = hist(l[[i]]$probability, breaks=hists[[1]]$breaks, col=col[i], add=TRUE)
		hists = c(hists,list(h))
		
		# Label accuracy
		if(labels && i%%2==0){
			if(!is.na(hists[[i-1]]$counts) && length(hists[[i-1]]$counts > 0)){
				lbl = hists[[i]]$counts/hists[[i-1]]$counts
				lbl[which(lbl %in% NaN | lbl %in% 0)] = ""
				lbl = sprintf("%.2f",lbl)
				lbl[which(lbl %in% "NA")] = ""
				#lbl = format(lbl, digits=1)
				text(hists[[i-1]]$mids,hists[[i-1]]$count,labels=lbl,pos=3,cex=1)
			}
			# Put label on top of correct bar
			#text(hists[[i]]$mids,hists[[i]]$count,labels=lbl,pos=3,cex=1)
			
			# Label yield
			#lbl = hists[[i-1]]$counts;
		}
	}
}

rmsd_overlay = function(rmsd, good, xlab="RMSD", ylab="Count", title="RMSD to Native Structure"){
	b = hist(rmsd, xlab=xlab, ylab=ylab, main=title, breaks=20, col="gray")
	g = hist(good, breaks=b$breaks, col="seagreen", add=TRUE)
	labels = g$counts/b$counts;
	labels = format(labels, digits=3)
	text(b$mids,b$count,labels=labels,pos=3,cex=1)
}

plot_fp = function(data,
		main="False Positive Analysis",
		since_solved=T,
		adjusted=F,
		correct=T,
		xlim=c(0,1),
		ylim=c(0,length(data$id))){
		
		data_c = correct(data, since_solved=since_solved, adjusted=adjusted) 
		data_i = data[which(!data$id %in% data_c$id),]
		
		all = hist(data$prob, breaks=40, plot=F)
		inc = hist(data_i$prob,plot=F,breaks=all$breaks)
		cor = hist(data_c$prob,plot=F,breaks=all$breaks)
		ratio = cor$counts/all$counts
		
		plot(all$breaks[-1],inc$counts,type="l",col="black", xlab=NA, ylab=NA)
		par(new=T)
		plot(all$breaks[-1],cor$counts,type="l",col="red", xlab=NA,ylab=NA,axes=F)
		par(new=T)
		plot(all$breaks[-1],ratio,type="l",col="blue",ylim=c(0,1),xlim=c(0,1),xlab=NA,ylab=NA,axes=F)
		axis(4)
		title(main=main,xlab="MCM Score",ylab="Frequency",sub="Line plot of histogram bins")
		legend(0,1, c("Incorrect","Correct","Accuracy"), lty=1, lwd=2, col=c("black","red","blue"))
		return(list(all$breaks[-1],ratio))
}

plot_b = function(data,
		main="False Positive Analysis",
		since_solved=F,
		adjusted=F,
		correct=T,
		xlim=c(0,1),
		ylim=c(0,length(data$id)),
		ratio_xlim=xlim,
		ratio_ylim=c(0,1),
		ratio_int=(ratio_xlim[2]-ratio_xlim[1])/1000,
		count_xlim=xlim,
		count_ylim=ylim,
		count_int=(count_xlim[2]-count_xlim[1])/1000){

	if(correct){
		data_c = correct(data, since_solved=since_solved, adjusted=adjusted)
		data_i = data[which(!data$id %in% data_c$id),]
	} else {
		data_c = NA
		data_i = data
	}
	
	yield = function(x){
		return(length(which(data$prob >= x)))
	}
	
	cor = function(i){
		return(length(which(data_c$prob >= i)))
	}
	cor_v = function(x){
		count = c()
		for(i in x){
			count=c(count,cor(i))# & data_c$prob < i+count_int)))#/length(ssa1$prob))
		}
		return(count)
	}
	
	inc = function(i){
		return(length(which(data_i$prob >= i))) # & data_i$prob < i+count_int)))#/length(ssa1$prob))
	}
	inc_v = function(x){
		count = c()
		for(i in x){
			count=c(count,inc(i))
		}
		return(count)
	}
	
	ratio = function(i){
		d = length(which(data$prob >= i))
		if(correct){
			n = length(which(data_c$prob >= i))
		} else {
			n = accuracy(i)*d 
		}
		if(d==0){
			return(1)
		}
		return(n/d)
	}
	ratio_v = function(x){
		count = c()
		for(i in x){
			count=c(count,ratio(i))
		}
		return(count)
	}
	
	accuracy = function(x) {
		return(0.2994 + -0.5461*x + 0.9615*(x^2))
	}

	#if(correct){
		r = uniroot(function(x){ratio(x)-0.5}, interval=c(0,1))
		str(r)
		y = yield(r[[1]])
		l = length(data$prob)
		print("Yield: "||y||"/"||l||":"||sprintf("%0.3f",y/l))
		
		x_ratio = seq(ratio_xlim[1],ratio_xlim[2],ratio_int)
		xs = x_ratio
		y_ratio = ratio_v(xs)
		if(correct){
			plot(xs, y_ratio, type="l", col="blue", xlim=ratio_xlim,ylim=ratio_ylim, axes=F, ylab=NA, xlab=NA)
			par(new=T)
			points(r[[1]],0.5,pch=23)
			axis(3, col="blue", col.ticks="blue")
			axis(4, col="blue", col.ticks="blue")
		}
		#mtext("Accuracy", side=4, line=3, cex.lab=1,las=2, col="blue")
		#mtext("MCM Score", side=3, line=3, cex.lab=1, col="blue")
		

	#}
	
	xs = seq(count_xlim[1],count_xlim[2],count_int)
	#if(correct){
		
		if(correct){
			par(new=TRUE)
			plot(xs, cor_v(xs), type="l", col="red", xlim=count_xlim, ylim=count_ylim, axes=F, ylab=NA, xlab=NA)
		} else {
			plot(xs, inc_v(xs)*accuracy(xs), type="l", col="red", xlim=count_xlim, ylim=count_ylim, axes=F, ylab=NA, xlab=NA)
		}
		
		par(new=TRUE)
	#}
	if(correct){
		plot(xs, inc_v(xs), type="l", col="black", xlim=count_xlim, ylim=count_ylim, axes=F, ylab=NA, xlab=NA)
	} else {
		plot(xs, inc_v(xs)*(1-accuracy(xs)), type="l", col="black", xlim=count_xlim, ylim=count_ylim, axes=F, ylab=NA, xlab=NA)
	}
	par(new=T)
	mean_length = mean(data$length)
	if(!is.na(mean_length)){
		mean_length=" AvgSeqLen:"||sprintf("%0.0f",mean_length)
	} else {
		mean_length=""
	}
	title(xlab="MCM Score", ylab="Count", main=main||mean_length)
	axis(1)
	if(correct){
		axis(2)
	} else {
		axis(4)
	}
	if(correct){
		return(list(x_ratio,y_ratio))
	}
	title(sub="50/50 at MCM: "||sprintf("%0.2f",r[[1]])||" Yield: "||y||"/"||l||": "||sprintf("%0.3f",y/l))
	
}
