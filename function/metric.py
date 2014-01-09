import MySQLdb
from predictor import get_not_id
from predictor import is_not_id
from predictor import Predictor
from predictor import ALL_TERM
import cPickle
import math

TINY_NUM = 0.0000001
MAX_UPLOAD_LIST = 50

class Metric():
	"""
	Interface for certain metrics like log ratio and mutual information
		between two predictors.  Keys are stored as tuples.
	"""
	
	def __init__(self, dict=None, default=None):
		"""
		@param dict: Required dictionary type to support persistence.
		@type x: dict.
		@param default: If specified, the default value is returned instead of
			throwing key errors.
		"""
		#kdrew: for some weird reason you cannot properly initialize a dict in the parameter declaration
		if None == dict:
			dict = {}

		self._dict = dict
		self._default = default

	def __iter__(self):
		"""
		@return: iterator over the set of string tuple keys (acc1,acc2) known by
			this metric.
		"""
		for key in self._dict.__iter__():
			yield self._unkey(key)

	def get_all_ids(self, no_nots=True):
		"""
		@deprecated: Use iter instead.
		@return: Set of string tuple keys (acc1,acc2)
		"""
		if no_nots:
			ret_keys = []
			for key in self._dict.keys():
				key = self._unkey(key)
				#print "key: ", key, "key[0], key[1]: ", key[0], key[1]
				if not is_not_id(key[0]) and not is_not_id(key[1]):
					yield key
		else:
			for key in  self:
				yield key

	#kdrew: I think acc2 is suppose to be default to None, for some reason was not None, changed 8/14/09
	def get_metric(self, acc1, acc2=None):
		"""
		Retrieves the metric value for a metric M(acc2|acc1).
			If acc2 is None, the base metric is returned. 
		@param acc1: The observation.
		@type acc1: basestring or Predictor
		@param acc2: Optional. The target of your conditional metric.
		@type acc2:   basestring or Predictor
		"""
		assert(acc1 != None)
		key = self._key(acc1, acc2)
		try:
			return self._dict[key]
		except KeyError:
			#print "DEFAULT:",self._default
			if self._default != None:
				return self._default
			else:
				raise

	#kdrew: put in so I can get all metrics with a given acc
	def get_all_metric(self, acc1):
		"""
		Retrieves a dictionary of all metric values for a metric M(Everything|acc1).
		@param acc1: The observation.
		@type acc1: basestring or Predictor
		"""
		assert(acc1 != None)
		ret_dict = Metric()
		for given, target in self:
			if given == acc1 and target !=None:	
				ret_dict.set_metric(acc1=target, acc2=None,value=self.get_metric(given,target))
		return ret_dict
		
	def set_metric(self, acc1, acc2, value):
		"""
		Set the value of a metric.
		@return: old value or None. 
		"""
		assert(acc1 != None)
		old = None
		key = self._key(acc1, acc2)
		if self.contains_metric(acc1, acc2):
			old = self._dict[key]
		self._dict[key] = value
		return old

	def contains_metric(self, acc1, acc2):
		"""
		@return: True/False, dictionary contains this key combo.
		"""
		key = self._key(acc1, acc2)
		return self._dict.__contains__(key)

	def upload_metric(self, conn, table_name, delete_table=False):
		cursor = conn.cursor(MySQLdb.cursors.DictCursor)

		if delete_table:
			delete_query = "DROP TABLE IF EXISTS " + table_name
			cursor.execute(delete_query)

		create_query = "CREATE TABLE IF NOT EXISTS " + table_name + " ( id INT(11) NOT NULL auto_increment, acc VARCHAR(55), acc2 VARCHAR(55), metric DOUBLE default 0, PRIMARY KEY (id), KEY (acc), KEY (acc2) )"
		print create_query
		cursor.execute(create_query)

		insert_query = """INSERT INTO """ + table_name + """ (acc, acc2, metric) VALUES (%s,%s,%s)"""
		upload_list = []
		for acc, acc2 in self.get_all_ids():
			k = [acc, acc2]
			k.append(self.get_metric(acc, acc2))
			upload_list.append(k)

			#kdrew: flush every so often to avoid sql packet size error
			if len(upload_list) > MAX_UPLOAD_LIST:
				cursor.executemany(insert_query, upload_list)
				upload_list = []

		cursor.executemany(insert_query, upload_list)

	def load_metric(self, conn, table_name, 
				    query="select m.acc, m.acc2, m.metric from %s as m "):
		"""
		Load in the metrics from a database.
		@param table_name: The table where the metric exists.
		@param query: The selection query can be overriden.  The table_name is
			inserted into this query with query % table_name
		"""
		query = query % table_name
		print query
		cursor = conn.cursor()
		cursor.execute(query)
		print "Fetching"
		for acc1, acc2, metric in cursor.fetchall():
			self.set_metric(acc1, acc2, metric)
			
	def printTables(self):
		"""
		@deprecated: We shouldn't be printing any of these huge tables.
		"""
		print self._dict
		
	def _ids(self, acc1, acc2):
		"""
		Handles converting between Predictor and String types.  Should not be
			overriden.  To override key management @see: _key
		"""
		if isinstance(acc1, Predictor):
			acc1 = acc1.get_id()
		if isinstance(acc2, Predictor):
			acc2 = acc2.get_id()
		return (acc1, acc2)

	def _key(self, acc1, acc2):
		"""
		Stores keys as strings for shelve persistence.
		"""
		acc1,acc2 = self._ids(acc1, acc2)
		return ",".join([str(a) for a in (acc1, acc2)])
	
	def _unkey(self, key):
		"""
		Splits key strings back into tuples.
		"""
		all = key.split(",")
		ret = []
		for a in all:
			if a == 'None':
				a = None
			ret.append(a)
		return tuple(ret)
	
#kdrew: temporary checking to see if storing all keys as string works
#	def _key(self, acc1, acc2):
#		"""
#		Provides an override mechanism for persistent metric dictionaries to
#			maintain keys.
#		"""
#		return self._ids(acc1,acc2)
#	
#	def _unkey(self, key):
#		"""
#		Provides an override mechanism for persistent metric dictionaries to
#			maintain keys.
#		"""
#		return tuple(key)

class MetricShelve(Metric):
	"""
	Metric class that utilizes a shelve for its dictionary to conserve 
		memory and handle persistence.
	"""
	def __init__(self, filename=None, default=None):
		"""
		@param filename: Filename to use for the shelve.  If none a tempfile is
			created and deleted upon @see: close()
		"""
		Metric.__init__(self, dict=None, default=default)
		import shelve, tempfile, os
		if filename:
			self._temp = None
			self._filename = filename
		else:
			self._temp = tempfile.NamedTemporaryFile()
			self._filename = self._temp.name
			#kdrew: close the temporary file, just needed the temporary name
			self._temp.close()
		self._dict = shelve.open(self._filename,"n")
		
	def load_metric(self, conn, table_name, 
				    query="select m.acc, m.acc2, m.metric from %s as m "):
		"""Limits query to only conditional molecular functions."""
		query = query+" join mygolite_062009.term t on m.acc=t.acc and t.term_type='molecular_function'" 
		Metric.load_metric(self, conn, table_name, query)
	
	def close(self):
		self._dict.close()
		if self._temp:
			import os
			os.remove(self._filename)
			
	def _key(self, acc1, acc2):
		"""
		Stores keys as strings for shelve persistence.
		"""
		acc1,acc2 = self._ids(acc1, acc2)
		return ",".join([str(a) for a in (acc1, acc2)])
	
	def _unkey(self, key):
		"""
		Splits key strings back into tuples.
		"""
		all = key.split(",")
		ret = []
		for a in all:
			if a == 'None':
				a = None
			ret.append(a)
		return tuple(ret)

class MetricDB(Metric):
	"""Metric that queries at runtime."""
	
	def __init__(self, connection, table_name, default=None):
		Metric.__init__(self, default=default)
		self.connection = connection
		self.table_name = table_name
		
	def _load_metric(self, acc, acc2):
		"""
		Performs a query to load the metric into the db.
		"""
		def st(acc):
			if acc is None:
				return "is NULL"
			else:
				return "='%s'" % str(acc)
		
		acc,acc2 = self._ids(acc1, acc2)	
		query = """select m.metric from """ + self.table_name + """ as m where acc %s and acc2 %s""" % (st(acc), st(acc2))
		#print query
		cursor = self.connection.cursor()
		cursor.execute(query)
		metric = cursor.fetchone()[0]
		cursor.close()
		self.set_metric(acc, acc2, metric)
	
	def get_metric(self, acc1, acc2=None):
		if not self.contains_metric(acc1, acc2):
			try:
				self._load_metric(acc1, acc2)
			except:
				pass
		return Metric.get_metric(self, acc1, acc2)

class Frequency(Metric):
	def get_metric(self, acc1, acc2=None):
		if None == acc2:
			try:
				if is_not_id(acc1):
					all_freq = self.get_metric(ALL_TERM) 
					return  all_freq - self.get_metric(get_not_id(acc1.get_id()))
				else:
					return Metric.get_metric(self,acc1, None)
			except KeyError:
				return 0

		else:
			#kdrew: F(notacc2 , notacc1)
			if is_not_id(acc1) and is_not_id(acc2):
				acc_acc2_freq = self.get_metric(get_not_id(acc1), get_not_id(acc2))
				acc_freq = self.get_metric(get_not_id(acc1))
				acc2_freq = self.get_metric(get_not_id(acc2))
				all_freq = self.get_metric(ALL_TERM) 
				return all_freq - acc_freq - acc2_freq + acc_acc2_freq

			#kdrew: F(notacc2 , acc1)
			elif is_not_id(acc2):
				acc_acc2_freq = self.get_metric(acc1, get_not_id(acc2))
				acc_freq = self.get_metric(acc1)
				return acc_freq - acc_acc2_freq

			#kdrew: F(acc2 , notacc1)
			elif is_not_id(acc1):
				acc_acc2_freq = self.get_metric(acc2, get_not_id(acc1))
				acc2_freq = self.get_metric(acc2)
				return acc2_freq - acc_acc2_freq

			#kdrew: F(acc1, acc2)
			else:
				try:
					return Metric.get_metric(self,acc1, acc2)

				except KeyError:
					#kdrew: if the one combination of predictors doesn't work reverse them and try again
					try:
						return Metric.get_metric(self,acc2, acc1)
					except KeyError:
						return 0

	def _increment_freq_metric(self, acc1, acc2=None):
		if self.contains_metric(acc1,acc2):
			self.set_metric(acc1, acc2, self.get_metric(acc1,acc2)+1)
		else:
			#print "new key: ", acc1, acc2
			self.set_metric(acc1, acc2,1)
		
	def compute_metric_cluster(self, cluster):
		#kdrew: for all acc in clstr entry
		for clstr_pred in cluster:
			#kdrew: add count to just acc 
			self._increment_freq_metric(clstr_pred.get_id())

			#kdrew: for all preds in clstr entry
			for clstr_pred2 in cluster:
				#kdrew: add to pred pred pairs
				#print "clstr_pred: ",clstr_pred.get_id(), " clstr_pred2: ", clstr_pred2.get_id()
				self._increment_freq_metric(clstr_pred.get_id(), clstr_pred2.get_id())

	#kdrew: this function computes frequency tables
	#kdrew: clstrs: dictionary of clusters where keys are mygo seq_ids and entries is a list of terms,
	def compute_metric(self, clstrs):
		#kdrew: for all sequence ids in dictionary of accs
		for clstr in clstrs:
			print "clstr_id: ",clstr
			self.compute_metric_cluster(clstrs[clstr])


class LogRatios(Metric):

	DEFAULT_MIN = math.log(TINY_NUM)

	def __init__(self, dict=None, default=None):
		Metric.__init__(self, dict=dict, default=LogRatios.DEFAULT_MIN)

	def compute_metric(self, prob_metric):
		"""
		Compute log ratios from the probability table.
		"""
		for acc,acc2 in prob_metric:

			notacc = get_not_id(acc)
			#kdrew: for all acc2's in prob_metric
			notacc2 = get_not_id(acc2)

			acc2_acc = prob_metric.get_metric(acc, acc2)
			acc2_notacc = prob_metric.get_metric(notacc, acc2)
			lr = acc2_acc / (acc2_notacc + TINY_NUM)
			#print "likelihood ratio: ", lr
			#kdrew: compute log ratio and add in TINY_NUM to avoid log(0)
			log_value = math.log(lr + TINY_NUM)
			self.set_metric(acc, acc2,log_value)



class LogRatiosShelve(LogRatios,MetricShelve):

	def __init__(self, filename=None):
		MetricShelve.__init__(self,
							  filename=filename,
							  default=LogRatios.DEFAULT_MIN)
		
class MutualInformation(Metric):

	def __init__(self,dict=None,default=0):
		Metric.__init__(self, dict=dict, default=default)

	def compute_metric(self, prob_metric):
		"""
		Computes mutual information from the probability tables
		@param prob_metric: The probability metrics to compute MI for. 
		@type prob_metric: Probability.
		"""
		for acc,acc2 in prob_metric:
			mi = self.compute_mi(acc, acc2, prob_metric)
			self.set_metric(acc, acc2, mi)

	def compute_mi(self, acc, acc2, prob_metric):
		acc_prob = prob_metric.get_metric(acc,None)
		acc2_prob = prob_metric.get_metric(acc2,None)
		not_acc_prob = prob_metric.get_metric(get_not_id(acc),None)
		not_acc2_prob = prob_metric.get_metric(get_not_id(acc2),None)

		acc_acc2_prob = prob_metric.get_metric(acc, acc2) * acc_prob
		not_acc_acc2_prob = prob_metric.get_metric(get_not_id(acc), acc2) * not_acc_prob
		acc_not_acc2_prob = prob_metric.get_metric(acc, get_not_id(acc2)) * acc_prob
		not_acc_not_acc2_prob = prob_metric.get_metric(get_not_id(acc), get_not_id(acc2)) * acc_prob

		tmp_MI = acc_acc2_prob * math.log((acc_acc2_prob / (acc_prob * acc2_prob + TINY_NUM)) + TINY_NUM)
		tmp_MI += not_acc_acc2_prob * math.log((not_acc_acc2_prob / (not_acc_prob * acc2_prob + TINY_NUM)) + TINY_NUM)

		#kdrew: generalize this so it can be a parameter (class level)
		#kdrew: original code only did P(acc2|acc) and P(acc2|not_acc)
		#kdrew: do not do "not" "not" because general terms wash out everything
		#tmp_MI += acc_not_acc2_prob * math.log((acc_not_acc2_prob/(acc_prob * not_acc2_prob+TINY_NUM))+TINY_NUM)
		#tmp_MI += not_acc_not_acc2_prob * math.log((not_acc_not_acc2_prob/(not_acc_prob * not_acc2_prob+TINY_NUM))+TINY_NUM)

		return tmp_MI

	def max(self, function, terms):
		max_term = None
		max_mi = -99

		for term in terms:
			#term.print_term()
			#print self.mutual_infos[(term.acc,function.acc)]
			if self.get_metric(function, term) >= max_mi:
				max_term = term
				max_mi = self.get_metric(function,term)
				#print "new max!"
		#if max_mi!=0:
		#	print "MI"
		return max_term

class MutualInformationShelve(MutualInformation,MetricShelve):

	def __init__(self, filename=None):
		MetricShelve.__init__(self,
							  filename=filename,
							  default=0)
	
#kdrew: conditional probability class
class Probability(Metric):
	def __init__(self, pc=0, dict=None):
		Metric.__init__(self,dict=dict)
		self.pseudo_count_value = pc
		self.hasPseudoCounts = False
		#print "pseudo_count:", self.pseudo_count_value


	#kdrew: order matters, ie. P(acc2|acc1)
	def get_metric(self, acc1, acc2):
		try:
			if acc1 == None:
				return TINY_NUM
			if acc2 == None:
				if is_not_id(acc1):
					return (1.0 - self.get_metric(get_not_id(acc1),None))
				else:
					return self._dict[self._key(acc1, None)]
			else:
				#kdrew: P(notacc2 | notacc1)
				if is_not_id(acc1) and is_not_id(acc2):
					rawacc2 = get_not_id(acc2)
					#kdrew: P(notacc2|notacc1)
					return 1.0 - self.get_metric(acc1, rawacc2)

				#kdrew: P(notacc2 | acc1)
				elif is_not_id(acc2):
					rawacc2 = get_not_id(acc2)
					#kdrew: P(notacc2|acc1)
					return 1.0 - self.get_metric(acc1, rawacc2)

				#kdrew: P(acc2 | notacc1)
				elif is_not_id(acc1):
					rawacc1 = get_not_id(acc1)
					p1 = self.get_metric(rawacc1,None)
					p2 = self.get_metric(acc2,None)
					p2G1 = self.get_metric(rawacc1, acc2)

					#kdrew: tests for small numerator and returns TINY_NUM if smaller
					#kdrew: fixes problem with P(all|notall) returning 1.0
					if (p2 - p2G1 * p1) <= TINY_NUM:
						return TINY_NUM
					else:
						#kdrew: P(acc2|notacc1)
						return (p2 - p2G1 * p1) / (1.0 - p1)

				else:
					return self._dict[self._key(acc1, acc2)]

		except KeyError:
			return TINY_NUM



	#kdrew: this function computes probability tables
	def compute_metric(self, freq_metric, pseudo_count=False):
		all_keys = freq_metric.get_all_ids()

		#kdrew: for all keys in freq_metric
		for key in all_keys:
			if isinstance(key, tuple):
				acc = key[0]
				acc2 = key[1]
			else:
				acc = key
				acc2 = None

			#self._dict[self._key((acc, acc2))] = self.compute_prob(freq_metric, acc, acc2)
			self.set_metric(acc, acc2, self.compute_prob(freq_metric, acc, acc2))

			if pseudo_count and None != acc2:
				self.compute_pseudo_count(freq_metric, acc, acc2)


	def compute_pseudo_count(self, freq_metric, acc, acc2):
		weighted_prob = self.get_metric(acc, acc2)
		weighted_prob = weighted_prob * freq_metric.get_metric(acc, acc2) 

		#kdrew: calculate background probability from frequency tables
		weighted_pc_prob = (freq_metric.get_metric(acc2) / (freq_metric.get_metric(ALL_TERM) + TINY_NUM)) * self.pseudo_count_value
		normalization = freq_metric.get_metric(acc, acc2) + self.pseudo_count_value + TINY_NUM

		self.set_metric(acc,  acc2,(weighted_prob + weighted_pc_prob) / normalization)

	def compute_prob(self, freq_metric, acc, acc2=None):
		if acc2 == None:
			#kdrew: P(acc)
			return freq_metric.get_metric(acc) / (freq_metric.get_metric(ALL_TERM) + TINY_NUM)
		else:
			#kdrew: P(acc2|acc)
			return freq_metric.get_metric(acc, acc2) / (freq_metric.get_metric(acc) + TINY_NUM)





