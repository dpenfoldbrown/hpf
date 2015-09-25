
library(ROCR)
source("/Users/kdrew/scripts/function_prediction/sql_setup.R")

data = ss_qry("SELECT mc.cluster_id AS cluster_id, SUBSTRING_INDEX(mc.experiment_sccs,\'.\',3) AS pred_sf, SUBSTRING_INDEX(ss.sccs,\'.\',3) as true_sf, mc.probability as top_prob,
		if(SUBSTRING_INDEX(mc.experiment_sccs,\'.\',3) = SUBSTRING_INDEX(ss.sccs,\'.\',3),1,0) as correct
		FROM since_solved.mcm_cluster mc, since_solved.since_solved ss, (SELECT cluster_id, MAX(probability) as max_prob from since_solved.mcm_cluster group by cluster_id) as max_table
		WHERE mc.cluster_id = ss.cluster_id and max_table.cluster_id = mc.cluster_id and mc.probability = max_table.max_prob
			AND ss.ascession_date >= DATE(\'2005-01-01\')
			AND mc.probability > 0
			AND mc.cluster_id not in (10224, 20062 ,26089 ,19291 ,11937 ,12920 ,20519 ,29906 ,18292 ,17200 ,18446 ,29715 ,29093 ,31091 ,29928 ,36059 ,11851 ,49033 ,51070 ,2241 ,36144 ,24456 ,8946 ,9534 ,43224 ,42776 ,17514 ,50872 ,11836 ,19018 ,42742 ,42346 ,31919 ,41777 ,43218 ,37759 ,30191 ,18886 ,30166 ,47095 ,13212 ,30173 ,45454 ,23176 ,44478 ,44991 ,45474 ,44978 ,13768 ,45933 ,7943)
		ORDER BY top_prob")


data2 = ss_qry("SELECT distinct cc.cluster_id AS cluster_id, SUBSTRING_INDEX(mcr.experiment_sccs,\'.\',3) AS pred_sf, SUBSTRING_INDEX(ss.sccs,\'.\',3) AS true_sf, mcr.probability as top_prob,
			if(SUBSTRING_INDEX(mcr.experiment_sccs,\'.\',3) = SUBSTRING_INDEX(ss.sccs,\'.\',3),1,0) as correct
		FROM since_solved.cdhit_clstr cc, since_solved.since_solved ss, since_solved.mcmData_redux mcr,
			(SELECT cc.cluster_id, MAX(mcr.probability) AS max_prob FROM since_solved.cdhit_clstr cc, since_solved.mcmData_redux as mcr WHERE mcr.sequence_key = cc.foreign_key group by cc.cluster_id) as max_table
		WHERE cc.cluster_id = ss.cluster_id and max_table.cluster_id = cc.cluster_id and mcr.probability = max_table.max_prob and mcr.sequence_key = cc.foreign_key
			AND ss.ascession_date >= DATE(\'2005-01-01\')
			AND mcr.probability > 0
			AND cc.cluster_id not in (10224, 20062 ,26089 ,19291 ,11937 ,12920 ,20519 ,29906 ,18292 ,17200 ,18446 ,29715 ,29093 ,31091 ,29928 ,36059 ,11851 ,49033 ,51070 ,2241 ,36144 ,24456 ,8946 ,9534 ,43224 ,42776 ,17514 ,50872 ,11836 ,19018 ,42742 ,42346 ,31919 ,41777 ,43218 ,37759 ,30191 ,18886 ,30166 ,47095 ,13212 ,30173 ,45454 ,23176 ,44478 ,44991 ,45474 ,44978 ,13768 ,45933 ,7943)
		ORDER BY top_prob")


pred <- prediction(data$top_prob, data$correct)
perf <- performance(pred,"prec", "rpp")

plot(perf, colorize=T, main="SCOP 1.67", lwd=3)
quartz()

pred2 <- prediction(data2$top_prob, data2$correct)
perf2 <- performance(pred2,"prec", "rpp")
plot(perf2, colorize=T, main="SCOP 1.75", lwd=3)

quartz()
plot(perf, colorize=T, lwd=2)
plot(perf2, colorize=T, lwd=2, add=T)

print("50% accuracy mcm value")
print(perf@alpha.values[[1]][max(which(perf@y.values[[1]] > .5))])

print("50% accuracy mcm value, SCOP 1.75")
print(perf2@alpha.values[[1]][max(which(perf2@y.values[[1]] > .5))])

print("50% accuracy mcm value, SCOP 1.75")
print(perf2@alpha.values[[1]][max(which(perf2@y.values[[1]] > .5))])

