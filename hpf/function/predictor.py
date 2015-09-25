
ALL_TERM = "all"
#kdrew: declaration of how to put negate a term's id
#kdrew: generally used for probability dictionaries
NOT = "not"

#kdrew: returns opposite acc 
#kdrew: nota.51.1 -> a.51.1 -> nota.51.1
def get_not_id(acc):
	"""
	@type pred: Predictor or basestring 
	@return: The negation of a term
	"""
	if None == acc:
		return None
	elif isinstance(acc, basestring):
		if is_not_id(acc):
			return acc.replace(NOT,"")
		else:
			return NOT + acc
	else:
		id = acc.get_id()
		assert(isinstance(id, basestring))
		return Predictor(get_not_id(id))
		
def is_not_id(pred):
	"""
	@type pred: Predictor or basestring 
	@return: True/False this term is a negation
	"""
	if None == pred:
		return False
	elif isinstance(pred, basestring):
		return pred.find("not") > -1
	else:
		id = pred.get_id()
		assert(isinstance(id, basestring))
		return is_not_id(id)

#kdrew: class which is generic for terms and superfamilies
class Predictor:

	def __init__(self, id=0):
		self.id = id

	def get_id(self):
		return self.id

	#kdrew: generally used for probability dictionaries
	def get_not_id(self):
		return NOT + self.get_id()

