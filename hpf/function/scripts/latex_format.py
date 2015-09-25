from __future__ import division
import string

class LatexFormatter():

	def denovo_byClass_table(self,a_n,a_n_corr,b_n,b_n_corr,c_n,c_n_corr,d_n,d_n_corr,o_n,o_n_corr):
		t_dict = {}
		t_dict['a_total'] = sum(a_n)
		t_dict['a_total_corr']= sum(a_n_corr)
		t_dict['a_high_conf']= a_n[-1:]
		t_dict['a_high_conf_corr'] = a_n_corr[-1:]
		t_dict['a_med_conf']= a_n[-2:-1] + t_dict['a_high_conf']
		t_dict['a_med_conf_corr'] = a_n_corr[-2:-1] + t_dict['a_high_conf_corr']

		t_dict['b_total'] = sum(b_n)
		t_dict['b_total_corr'] = sum(b_n_corr)
		t_dict['b_high_conf'] = b_n[-1:]
		t_dict['b_high_conf_corr'] = b_n_corr[-1:]
		t_dict['b_med_conf']= b_n[-2:-1] + t_dict['b_high_conf']
		t_dict['b_med_conf_corr'] = b_n_corr[-2:-1] + t_dict['b_high_conf_corr']

		t_dict['c_high_conf'] = c_n[-1:]
		t_dict['c_high_conf_corr'] = c_n_corr[-1:]
		t_dict['c_total'] = sum(c_n)
		t_dict['c_total_corr'] = sum(c_n_corr)
		t_dict['c_med_conf']= c_n[-2:-1] + t_dict['c_high_conf']
		t_dict['c_med_conf_corr'] = c_n_corr[-2:-1] + t_dict['c_high_conf_corr']

		t_dict['d_total'] = sum(d_n)
		t_dict['d_total_corr'] = sum(d_n_corr)
		t_dict['d_high_conf'] = d_n[-1:]
		t_dict['d_high_conf_corr'] = d_n_corr[-1:]
		t_dict['d_med_conf']= d_n[-2:-1] + t_dict['d_high_conf']
		t_dict['d_med_conf_corr'] = d_n_corr[-2:-1] + t_dict['d_high_conf_corr']

		t_dict['o_total'] = sum(o_n)
		t_dict['o_total_corr'] = sum(o_n_corr)
		t_dict['o_high_conf'] = o_n[-1:]
		t_dict['o_high_conf_corr'] = o_n_corr[-1:]
		t_dict['o_med_conf']= o_n[-2:-1] + t_dict['o_high_conf']
		t_dict['o_med_conf_corr'] = o_n_corr[-2:-1] + t_dict['o_high_conf_corr']

		return self.denovo_byClass_correct_table_format(t_dict)

	def denovo_byClass_correct_table_format(self,t_dict):
		total = t_dict['a_total']+t_dict['b_total']+ t_dict['c_total']+ t_dict['d_total']+t_dict['o_total']
		total_corr = t_dict['a_total_corr']+ t_dict['b_total_corr']+ t_dict['c_total_corr']+ t_dict['d_total_corr']+t_dict['o_total_corr']
		total_percent_corr = total_corr/total
		high_conf = t_dict['a_high_conf']+ t_dict['b_high_conf'] + t_dict['c_high_conf'] + t_dict['d_high_conf']+t_dict['o_high_conf']
		high_conf_corr = t_dict['a_high_conf_corr'] + t_dict['b_high_conf_corr'] + t_dict['c_high_conf_corr'] + t_dict['d_high_conf_corr'] +t_dict['o_high_conf_corr']
		high_conf_percent_corr = high_conf_corr/high_conf
		high_conf_percent = high_conf/total
		high_conf_percent_yield = high_conf/total

		med_conf = t_dict['a_med_conf']+ t_dict['b_med_conf'] + t_dict['c_med_conf'] + t_dict['d_med_conf']+t_dict['o_med_conf']
		med_conf_corr = t_dict['a_med_conf_corr'] + t_dict['b_med_conf_corr'] + t_dict['c_med_conf_corr'] + t_dict['d_med_conf_corr'] +t_dict['o_med_conf_corr']
		med_conf_percent_corr = med_conf_corr/med_conf
		med_conf_percent = med_conf/total
		med_conf_percent_yield = med_conf/total


		a_percent_corr = t_dict['a_total_corr']/t_dict['a_total']
		a_high_conf_percent_corr = t_dict['a_high_conf_corr']/t_dict['a_high_conf']
		a_percent = t_dict['a_total']/total
		a_high_conf_percent= t_dict['a_high_conf']/t_dict['a_total']
		a_high_conf_percent_yield = t_dict['a_high_conf']/total
		a_med_conf_percent_corr = t_dict['a_med_conf_corr']/t_dict['a_med_conf']
		a_med_conf_percent= t_dict['a_med_conf']/t_dict['a_total']
		a_med_conf_percent_yield = t_dict['a_med_conf']/total

		b_percent_corr = t_dict['b_total_corr']/t_dict['b_total']
		b_high_conf_percent_corr = t_dict['b_high_conf_corr']/t_dict['b_high_conf']
		b_percent = t_dict['b_total']/total
		b_high_conf_percent= t_dict['b_high_conf']/t_dict['b_total']
		b_high_conf_percent_yield = t_dict['b_high_conf']/total
		b_med_conf_percent_corr = t_dict['b_med_conf_corr']/t_dict['b_med_conf']
		b_med_conf_percent= t_dict['b_med_conf']/t_dict['b_total']
		b_med_conf_percent_yield = t_dict['b_med_conf']/total

		c_percent_corr = t_dict['c_total_corr']/t_dict['c_total']
		c_high_conf_percent_corr = t_dict['c_high_conf_corr']/t_dict['c_high_conf']
		c_percent = t_dict['c_total']/total
		c_high_conf_percent= t_dict['c_high_conf']/t_dict['c_total']
		c_high_conf_percent_yield = t_dict['c_high_conf']/total
		c_med_conf_percent_corr = t_dict['c_med_conf_corr']/t_dict['c_med_conf']
		c_med_conf_percent= t_dict['c_med_conf']/t_dict['c_total']
		c_med_conf_percent_yield = t_dict['c_med_conf']/total

		d_percent_corr = t_dict['d_total_corr']/t_dict['d_total']
		d_high_conf_percent_corr = t_dict['d_high_conf_corr']/t_dict['d_high_conf']
		d_percent = t_dict['d_total']/total
		d_high_conf_percent= t_dict['d_high_conf']/t_dict['d_total']
		d_high_conf_percent_yield = t_dict['d_high_conf']/total
		d_med_conf_percent_corr = t_dict['d_med_conf_corr']/t_dict['d_med_conf']
		d_med_conf_percent= t_dict['d_med_conf']/t_dict['d_total']
		d_med_conf_percent_yield = t_dict['d_med_conf']/total

		o_percent_corr = t_dict['o_total_corr']/t_dict['o_total']
		o_high_conf_percent_corr = t_dict['o_high_conf_corr']/t_dict['o_high_conf']
		o_percent = t_dict['o_total']/total
		o_high_conf_percent= t_dict['o_high_conf']/t_dict['o_total']
		o_high_conf_percent_yield = t_dict['o_high_conf']/total
		o_med_conf_percent_corr = t_dict['o_med_conf_corr']/t_dict['o_med_conf']
		o_med_conf_percent= t_dict['o_med_conf']/t_dict['o_total']
		o_med_conf_percent_yield = t_dict['o_med_conf']/total

		#kdrew: CLASS & class total (percent) & correct (percent) & med_conf (percent) & med_conf correct (percent) & med_conf yield percentage & high_conf (percent) & high_conf correct (percent) & high_conf yield percentage
		#kdrew: ALL & total & total correct (percent) & med_conf total (percent) & med_conf correct (percent) & & high_conf total (percent) & high_conf correct (percent) &

		table_str = """ 
                         A & %s (%s\\%%) &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & %s\\%% & %s (%s\\%%) & %s (%s\\%%) & %s\\%% \\\\
                         B & %s (%s\\%%) &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & %s\\%% & %s (%s\\%%) & %s (%s\\%%) & %s\\%% \\\\
                         C & %s (%s\\%%) &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & %s\\%% & %s (%s\\%%) & %s (%s\\%%) & %s\\%% \\\\
                         D & %s (%s\\%%) &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & %s\\%% & %s (%s\\%%) & %s (%s\\%%) & %s\\%% \\\\
                         Other & %s (%s\\%%) &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & %s\\%% & %s (%s\\%%) & %s (%s\\%%) & %s\\%% \\\\
			 \\hline
                         All & %s &  %s (%s\\%%) & %s (%s\\%%) & %s (%s\\%%) & & %s (%s\\%%) & %s (%s\\%%) & \\\\
		""" 

		return table_str % (int(t_dict['a_total']), round(a_percent,3)*100, int(t_dict['a_total_corr']), round(a_percent_corr,3)*100, 
						int(t_dict['a_med_conf']), round(a_med_conf_percent,3)*100, int(t_dict['a_med_conf_corr']), round(a_med_conf_percent_corr,3)*100, round(a_med_conf_percent_yield,3)*100, 
						int(t_dict['a_high_conf']), round(a_high_conf_percent,3)*100, int(t_dict['a_high_conf_corr']), round(a_high_conf_percent_corr,3)*100, round(a_high_conf_percent_yield,3)*100, 
					int(t_dict['b_total']), round(b_percent,3)*100, int(t_dict['b_total_corr']), round(b_percent_corr,3)*100, 
						int(t_dict['b_med_conf']), round(b_med_conf_percent,3)*100,int(t_dict['b_med_conf_corr']), round(b_med_conf_percent_corr,3)*100, round(b_med_conf_percent_yield,3)*100, 
						int(t_dict['b_high_conf']), round(b_high_conf_percent,3)*100,int(t_dict['b_high_conf_corr']), round(b_high_conf_percent_corr,3)*100, round(b_high_conf_percent_yield,3)*100, 
					int(t_dict['c_total']), round(c_percent,3)*100, int(t_dict['c_total_corr']), round(c_percent_corr,3)*100, 
						int(t_dict['c_med_conf']), round(c_med_conf_percent,3)*100,int(t_dict['c_med_conf_corr']), round(c_med_conf_percent_corr,3)*100, round(c_med_conf_percent_yield,3)*100, 
						int(t_dict['c_high_conf']), round(c_high_conf_percent,3)*100,int(t_dict['c_high_conf_corr']), round(c_high_conf_percent_corr,3)*100, round(c_high_conf_percent_yield,3)*100, 
					int(t_dict['d_total']), round(d_percent,3)*100, int(t_dict['d_total_corr']), round(d_percent_corr,3)*100, 
						int(t_dict['d_med_conf']), round(d_med_conf_percent,3)*100,int(t_dict['d_med_conf_corr']), round(d_med_conf_percent_corr,3)*100, round(d_med_conf_percent_yield,3)*100, 
						int(t_dict['d_high_conf']), round(d_high_conf_percent,3)*100,int(t_dict['d_high_conf_corr']), round(d_high_conf_percent_corr,3)*100, round(d_high_conf_percent_yield,3)*100, 
					int(t_dict['o_total']), round(o_percent,3)*100, int(t_dict['o_total_corr']), round(o_percent_corr,3)*100, 
						int(t_dict['o_med_conf']), round(o_med_conf_percent,3)*100,int(t_dict['o_med_conf_corr']), round(o_med_conf_percent_corr,3)*100, round(o_med_conf_percent_yield,3)*100,
						int(t_dict['o_high_conf']), round(o_high_conf_percent,3)*100,int(t_dict['o_high_conf_corr']), round(o_high_conf_percent_corr,3)*100, round(o_high_conf_percent_yield,3)*100,
					int(total), int(total_corr), round(total_percent_corr,3)*100, 
						int(med_conf), round(med_conf_percent,3)*100, int(med_conf_corr), round(med_conf_percent_corr,3)*100,
						int(high_conf), round(high_conf_percent,3)*100, int(high_conf_corr), round(high_conf_percent_corr,3)*100)


	def denovo_byClass_table_format_old(self,t_dict):

		total = t_dict['a_total']+t_dict['b_total']+ t_dict['c_total']+ t_dict['d_total']
		high_conf = t_dict['a_high_conf']+ t_dict['b_high_conf'] + t_dict['c_high_conf'] + t_dict['d_high_conf']
		a_high_conf_percent= t_dict['a_high_conf']/t_dict['a_total']
		b_high_conf_percent= t_dict['b_high_conf']/t_dict['b_total']
		c_high_conf_percent= t_dict['c_high_conf']/t_dict['c_total']
		d_high_conf_percent= t_dict['d_high_conf']/t_dict['d_total']

		a_high_conf_percent_yield = t_dict['a_high_conf']/total
		b_high_conf_percent_yield = t_dict['b_high_conf']/total
		c_high_conf_percent_yield = t_dict['c_high_conf']/total
		d_high_conf_percent_yield = t_dict['d_high_conf']/total
		high_conf_percent_yield = high_conf/total

		a_percent = t_dict['a_total']/total
		b_percent = t_dict['b_total']/total
		c_percent = t_dict['c_total']/total
		d_percent = t_dict['d_total']/total


		table_str = """ 
                         A & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\
                         B & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\
                         C & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\
                         D & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\
			 \\hline
                         All & %s &  %s (%s\\%%)  & \\\\
		""" 

		return table_str % (int(t_dict['a_total']), round(a_percent,3)*100, int(t_dict['a_high_conf']), round(a_high_conf_percent,3)*100, round(a_high_conf_percent_yield,3)*100, int(t_dict['b_total']), round(b_percent,3)*100, int(t_dict['b_high_conf']), round(b_high_conf_percent,3)*100, round(b_high_conf_percent_yield,3)*100, int(t_dict['c_total']), round(c_percent,3)*100, int(t_dict['c_high_conf']), round(c_high_conf_percent,3)*100, round(c_high_conf_percent_yield,3)*100, int(t_dict['d_total']), round(d_percent,3)*100, int(t_dict['d_high_conf']), round(d_high_conf_percent,3)*100, round(d_high_conf_percent_yield,3)*100, int(total),int(high_conf), round(high_conf_percent_yield,3)*100)


	def denovo_byClass_table_format(self,t_dict, classes=('a','b','c','d'), others=('e','f','g','i','h','j','k')):

		total = 0
		high_conf = 0
		for key in t_dict.keys():
			if key.find("total") > -1:
				total += t_dict[key]
			if key.find("high_conf") > -1:
				high_conf += t_dict[key]

		table_str = ""
		for scop_class in classes:
			high_conf_percent = t_dict[scop_class+'_high_conf']/t_dict[scop_class+'_total']
			high_conf_percent_yield = t_dict[scop_class+'_high_conf']/total
			percent = t_dict[scop_class+'_total']/total

			table_str += "%s & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\ \n" % (string.capitalize(scop_class), int(t_dict[scop_class+'_total']), round(percent,3)*100, int(t_dict[scop_class+'_high_conf']), round(high_conf_percent,3)*100, round(high_conf_percent_yield,3)*100)

		other_hc = 0
		other_total = 0
		for scop_class in others:
			try:
				other_hc += t_dict[scop_class+'_high_conf']
				other_total += t_dict[scop_class+'_total']
			except KeyError:
				continue
		
		table_str += "Other & %s (%s\\%%)  &  %s (%s\\%%)  &  %s\\%% \\\\ \n" % (int(other_total), round(other_total/total,3)*100, int(other_hc), round(other_hc/other_total,3)*100, round(other_hc/total,3)*100)
		table_str += "\\hline \n"
		high_conf_percent_yield = high_conf/total
		table_str += "All & %s &  %s (%s\\%%)  & \\\\" % (int(total),int(high_conf), round(high_conf_percent_yield,3)*100)



		return table_str


