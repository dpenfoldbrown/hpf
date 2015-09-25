
import unittest
import os
from scipy import special
from math import exp

#kdrew: taken from http://www.velocityreviews.com/forums/t352824-hypergeometric-distribution.html, Robert Kern
def logchoose(n, k):
	lgn1 = special.gammaln(n+1)
	lgk1 = special.gammaln(k+1)
	lgnk1 = special.gammaln(n-k+1)
	return lgn1 - (lgnk1 + lgk1)

def gauss_hypergeom(x, r, b, n):
	return exp(logchoose(r, x) + logchoose(b, n-x) - logchoose(r+b, n))


class SimpleTestCase(unittest.TestCase):
	def setUp(self):
		white_drawn = 4
		white_total = 5
		black_total = 45
		sampled = 10

		self.four_white = gauss_hypergeom(white_drawn, white_total, black_total, sampled)

		white_drawn = 5
		self.five_white = gauss_hypergeom(white_drawn, white_total, black_total, sampled)

	def testHG(self):
		assert self.four_white == 0.0039645830580152617, "incorrect hypergeometric distribution, four white"
		assert self.five_white == 0.00011893749174045687, "incorrect hypergeometric distribution, five white"


if __name__ == "__main__":
            unittest.main()
    


